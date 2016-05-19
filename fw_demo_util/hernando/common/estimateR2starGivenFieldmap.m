% Function: estimateR2starGivenFieldmap
%
% Description: estimate R2* map, given the fieldmap
% 
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = string containing the name of the species
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%% 
%  - fm: the estimated B0 field map
%
% Returns: 
%  - r2starmap: the estimated R2* map
%  - residual: fit error residual
%
% Author: Diego Hernando
% Date created: 2009
% Date last modified: August 18, 2011


function [r2starmap,residual] = estimateR2starGivenFieldmap ( imDataParams, algoParams, fm )


try
  precessionIsClockwise = imDataParams.PrecessionIsClockwise;
catch
  precessionIsClockwise = 1;
end

  
% If precession is clockwise (positive fat frequency) simply conjugate data
if precessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end


range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;
t = reshape(t,[],1);

sx = size(images,1);
sy = size(images,2);
C = size(images,4);
N = size(images,5);
num_acqs = size(images,6);

images = reshape(permute(images,[1 2 5 4 6 3]),[sx sy N C*num_acqs]);

% Undo effect of field map
images2 = zeros(size(images));
for kt=1:N
    for kc=1:C*num_acqs
      images2(:,:,kt,kc) = squeeze(images(:,:,kt,kc)).*exp(-j*2*pi*fm*t(kt));
    end
end

% Compute residual as a function of r2
r2s = linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
Phi = getPhiMatrixMultipeak(deltaF,relAmps,t);
P = [];
for k=1:NUM_R2STARS
  Psi = diag(exp(-r2s(k)*t));
  P = [P;(eye(N)-Psi*Phi*pinv(Psi*Phi))];
end

% Compute residual for all voxels and all field values
% Note: the residual is computed in a vectorized way, for increased speed
residual = zeros(sx,sy,NUM_R2STARS);

% Go line-by-line in the image to avoid using too much memory, while
% still reducing the loops significantly
% $$$ disp('Calculating residual...')
for kx=1:sx
  temp = reshape(squeeze(images2(kx,:,:,:)),[sy N*C*num_acqs]).';
  temp = reshape(temp,[N sy*C*num_acqs]);
  temp2 = reshape(sum(abs(reshape(P*temp,[N C*num_acqs*NUM_R2STARS*sy])).^2,1),[NUM_R2STARS C*num_acqs*sy]).';
  temp2 = sum(reshape(temp2,[C*num_acqs NUM_R2STARS*sy]),1);
  residual(kx,:,:) = reshape(temp2,[sy NUM_R2STARS]);
end
%  residual = shiftdim(residual,2);
% $$$ disp('done computing residual.');

[minres,iminres] = min(residual,[],3);
r2starmap = r2s(iminres);


