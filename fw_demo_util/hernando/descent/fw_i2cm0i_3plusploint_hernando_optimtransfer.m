% Function:  fw_i2cm0i_3plusploint_hernando_optimtransfer
%
% Description: Estimation of fieldmap, water image and fat image from
% a Dixon-type aquisition, based on the method proposed by Wuh and
% Fessler (ISMRM 2008), which uses Optimization Transfer to descend
% on the cost function
%
% (Huh W, Fessler JA, Samsonov AA. Water-fat decomposition with regularized field map. 
% In: Proceedings of the 16th Annual Meeting of ISMRM, Toronto, Canada, 2008. p. 1382.)
%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.fieldmap: initial fieldmap (in Hz) 
%%   - algoParams.r2starmap: initial R2* map (in 1/s) 
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%%     
%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%%   - algoParams.OT_ITERS = 40; % Number of optimization transfer iterations
%%   - algoParams.lambdamap ; % Spatially-varying regularization parameter
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
%
% Author: Diego Hernando
% Date created: June 11, 2009
% Date last modified: November 10, 2011

function outParams = fw_i2cm0i_3plusploint_hernando_optimtransfer( imDataParams, algoParams )

DEBUG_MODE = 0;

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% Get recon parameters and images
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
range_fm = algoParams.range_fm;
t = reshape(imDataParams.TE,[],1);
[sx,sy,sz,C,N,num_acqs] = size(imDataParams.images);
images = permute(imDataParams.images,[1 2 5 4 6 3]);

fm0 = algoParams.fieldmap;
r2starmap = algoParams.r2starmap;
NUM_ITERS = algoParams.OT_ITERS;
lambdamap = algoParams.lambdamap;

% Need to create my finite-difference matrix
if DEBUG_MODE==1
  tic
end
Dx = spalloc(sx*sy,sx*sy,2*sx*sy);
Dy = spalloc(sx*sy,sx*sy,2*sx*sy);
Dxy1 = spalloc(sx*sy,sx*sy,2*sx*sy);
Dxy2 = spalloc(sx*sy,sx*sy,2*sx*sy);
lmapvec = reshape(lambdamap,[],1);
for k=1:sx*sy
  if rem(k-1,sx) ~= 0
    Dx(k,k) = -min(lmapvec(k), lmapvec(k-1));
    Dx(k,k-1) = min(lmapvec(k), lmapvec(k-1));
  end
  if k-sx >= 1
    Dy(k,k) = -min(lmapvec(k), lmapvec(k-sx));
    Dy(k,k-sx) = min(lmapvec(k), lmapvec(k-sx));
  end
  if k-sx >= 1 & rem(k-1,sx) ~= 0
    Dxy1(k,k) = -min(lmapvec(k), lmapvec(k-sx-1));
    Dxy1(k,k-sx-1) = min(lmapvec(k), lmapvec(k-sx-1));
  end
  if k-sx >= 1 & rem(k+1,sx) ~= 0
    Dxy2(k,k) = -min(lmapvec(k), lmapvec(k-sx+1));
    Dxy2(k,k-sx+1) = min(lmapvec(k), lmapvec(k-sx+1));
  end
end
Dlambda = sqrt(2)*[Dx; Dy; Dxy1 ; Dxy2];
if DEBUG_MODE==1
  toc
end

% Need to get a upper bound on the Hessian of the data term 
% This should be easy, because data term is voxel-independent
Phi = getPhiMatrixMultipeak( deltaF,relAmps,t );
Gamma = Phi*inv(Phi'*Phi)*Phi';
[TM,TN] = meshgrid(reshape(t,[],1));
maxSignal = sum(sum(max(abs(images).^2,[],3),4),5);

if norm(r2starmap(:))==0
  B2 = 4*pi^2*sum(sum(abs(Gamma).*(TN-TM).^2));
  d2bound = B2*maxSignal;
else
  B2 = zeros(sx,sy);
  curGamma = zeros(sx,sy,N,N,num_acqs);
  for kx=1:sx
    for ky=1:sy
      for ka=1:num_acqs
        curPhi = Phi.*exp(-r2starmap(kx,ky,ka)*repmat(t(:),[1 2]));
        curGamma(kx,ky,:,:,ka) = curPhi*inv(curPhi'*curPhi)*curPhi';
        B2(kx,ky) = B2(kx,ky) + 4*pi^2*sum(sum(abs(squeeze(curGamma(kx,ky,:,:,ka))).*(TN-TM).^2));
      end
    end
  end
  d2bound = B2.*maxSignal;
end
  
D2data = spdiags(reshape(d2bound,[],1),0,sx*sy,sx*sy);
curHessTotal = Dlambda'*Dlambda + D2data;

fm = fm0;
for kit = 1:NUM_ITERS
  
  % Form the linear system... we need the gradient here
  d1 = zeros(sx,sy);
  for ka=1:num_acqs
    for kc=1:C
      for kn=1:N
        for km=1:N
          if norm(r2starmap(:))==0
            d1 = d1 - j*conj(images(:,:,kn,kc,ka)).*images(:,:,km,kc,ka)*Gamma(kn,km).*exp(j*2*pi*fm*(t(kn)-t(km)))*2*pi*(t(kn)-t(km));
          else
            d1 = d1 - j*conj(images(:,:,kn,kc,ka)).*images(:,:,km,kc,ka).*curGamma(:,:,kn,km,ka).*exp(j*2*pi*fm*(t(kn)-t(km)))*2*pi*(t(kn)-t(km));
          end
        end
      end
    end
  end
  d1 = real(d1); %Remove numerical errors that introduce imaginary component
  
  
  % Solve it
  curGradTotal = reshape(d1,[],1) + Dlambda'*Dlambda*fm(:);
  [fmstep,flag] = pcg(curHessTotal,-curGradTotal,1e-6,400);
  
  % Update field map
  fm = fm + reshape(fmstep,sx,sy);
  
  if DEBUG_MODE == 1
    curroughness = norm(Dlambda*fm(:))
    imagesc(reshape(fm,sx,sy),[-150 150]);drawnow
  end
end

% Now that we have the field map, estimate the water/fat images
for ka=1:num_acqs
  curParams = imDataParams;
  curParams.images = imDataParams.images(:,:,:,:,:,ka);
  amps = decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
  w(:,:,:,ka) = squeeze(amps(:,:,1,:));
  f(:,:,:,ka) = squeeze(amps(:,:,2,:));
end


% Put results in outParams structure
try
  outParams.species(1).name = algoParams.species(1).name;
  outParams.species(2).name = algoParams.species(2).name;
catch
  outParams.species(1).name = 'water';
  outParams.species(2).name = 'fat';
end  
  
outParams.species(1).amps = w;
outParams.species(2).amps = f;
outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;

