%% Function: graphCutIterations
%%
%% Description: graph cut iterations for field map estimation
%% 
%% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
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
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%   - algoParams.range_fm = [-400 400]; % Range of field map values
%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%%   - algoParams.lambda = 0.05; % Regularization parameter
%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
%%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 
%%
%%   - residual: the fit residual, of size NUM_FMS X sx X sy 
%%   - lmap: spatially varying regularization parameter
%%   - cur_ind: initial indices for field map (the indices determine the field map)
%%
%% Returns: 
%%  - fm: field map
%%
%% Author: Diego Hernando
%% Date created: June 3, 2009
%% Date last modified: August 18, 2011

function fm = graphCutIterations(imDataParams,algoParams,residual,lmap,cur_ind )

DEBUG = 0;
SMOOTH_NOSIGNAL = 1; % Whether to "homogenize" the lambdamap after
                     % some iterations, to get a smoother fieldmap in
                     % low-signal regions

STARTBIG = 1;

dkg = 15; % After dkg iterations, we may switch to a more homogeneous
          % regularization, to achieve more smoothness in noise-only
          % regions

DISPLAY_ITER = 0;


% Initialize some auxiliary variables
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:))*(imDataParams.FieldStrength)];
lambda = algoParams.lambda;
dt = imDataParams.TE(2)-imDataParams.TE(1);
period = 1/dt;
[sx,sy,N,C,num_acqs] = size(imDataParams.images);
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
resoffset = [0:(sx*sy-1)]'*algoParams.NUM_FMS;
[masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima( residual, 0.06 );
numLocalMin = size(resLocalMinima,1);
stepoffset = [0:(sx*sy-1)]'*numLocalMin;
ercurrent = 1e10;
clear cur_ind2
% $$$ [A] = createExpansionGraphVARPRO_fast( residual, dfm, ones(sx,sy), algoParams.size_clique,ones(sx,sy), ones(sx,sy));


% Main graph-cut loop. At each iteration, a max-flow/min-cut problem is solved
fm = zeros(sx,sy);
for kg=1:algoParams.NUM_ITERS

  fmiters(:,:,kg) = fm;
  
  if kg == 1 & STARTBIG==1
    lambdamap = lambda*lmap;
    ercurrent = 1e10;
    prob_bigJump = 1;
  elseif (kg==dkg & SMOOTH_NOSIGNAL==1) | STARTBIG == 0
    lambdamap = lambda*lmap;
    ercurrent = 1e10;
    prob_bigJump = 0.5;
  end
  
  % Get the sign of the current jump
  cur_sign = (-1)^(kg);
  
  % Configure the current jump move
  if rand < prob_bigJump % If we're making a "jumpmin" move.
    cur_ind2(1,:,:) = cur_ind;
    repCurInd = repmat(cur_ind2,[numLocalMin,1,1]);
    if cur_sign>0
      stepLocator = (repCurInd+20/dfm>=resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1))+1 ;
      validStep = masksignal>0 & stepLocator<=numMinimaPerVoxel;
    else
      stepLocator = (repCurInd-20/dfm>resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1));
      validStep = masksignal>0 & stepLocator>=1;
    end
    nextValue = zeros(sx,sy);
    nextValue(validStep) = resLocalMinima(stepoffset(validStep) + stepLocator(validStep));
    cur_step = zeros(sx,sy);
    cur_step(validStep) = nextValue(validStep) - cur_ind(validStep); 
    
    if rand < 0.5
      nosignal_jump = cur_sign*round(abs(deltaF(2))/dfm);
    else
      nosignal_jump = cur_sign*abs(round((period - abs(deltaF(2)))/dfm));  
    end
    cur_step(~validStep) = nosignal_jump;
    
  else % If we're making a standard jump move.
    all_jump = cur_sign*ceil(abs(randn*3));
    cur_step = all_jump*ones(sx,sy);
    nextValue = cur_ind + cur_step;
    
    % DH* 100309: fix errors where fm goes beyond range
    if cur_sign>0
      cur_step(nextValue(:)>length(fms)) = length(fms) - cur_ind(nextValue(:)>length(fms));
    else
      cur_step(nextValue(:)<1) = 1 - cur_ind(nextValue(:)<1);
    end
    
    if DEBUG == 1
      all_jump
    end
  end
  
  
  % Create the huge adjacency matrix
  if norm(lambdamap,'fro') > 0
    [A] = createExpansionGraphVARPRO_fast( residual, dfm, lambdamap, algoParams.size_clique,cur_ind, cur_step);
  else 
    [A] = createExpansionGraphVARPRO_fast( residual, dfm, lambdamap, algoParams.size_clique,cur_ind, cur_step);
  end
  A(A<0)=0;
  
  % Solve the max-flow/min-cut problem (max_flow is a function of the matlabBGL library)
  [flowvalTS,cut_TS,RTS,FTS] = max_flow(A',size(A,1),1);
    
  % Take the output of max_flow and update the fieldmap estimate with
  % the best neighbor (among a set of 2^Q neighbors, where Q is the
  % number of voxels)
  cut1 = (cut_TS==-1); 
  cut1b = 0*cut1 + 1;
  cut1b(end) = 0;
  if DEBUG == 1
    compareGC = [sum(sum(A(cut1b==1, cut1b==0)))  ,  sum(sum(A(cut1==1, cut1==0)))]
  end
  
  if sum(sum(A(cut1b==1, cut1b==0))) <= sum(sum(A(cut1==1, cut1==0))) 
    cur_indST = cur_ind;

    if DEBUG == 1
      disp('Not taken');
    end
  else
    cut = reshape(cut1(2:end-1)==0,sx,sy);
    cur_indST = cur_ind + cur_step.*(cut);
    erST = sum(sum(residual(cur_indST(:)+resoffset(:)))) + dfm^2*lambdamap(1,1)*sum(sum(abs(diff(cur_indST,2,1)).^2)) + dfm^2*lambdamap(1,1)*sum(sum(abs(diff(cur_indST,2,2)).^2));
    if DEBUG == 1
      disp('Taken');
% $$$   pixelschanged = sum(sum(cut));
    end
  end
  
  % Update the field map (described as a map of indices in this function)
  prev_ind = cur_ind; 
  cur_ind = cur_indST;

% $$$   size(fms)
% $$$   mymax = max(cur_ind(:))
  
  cur_ind(cur_ind<1) = 1;
  cur_ind(cur_ind>length(fms)) = length(fms);

  fm = fms(cur_ind);
  
  if DISPLAY_ITER == 1
    % If we want to see the fieldmap at each step
    imagesc(fm,[-600 600]);axis off equal tight,colormap gray;colorbar;
    title(['Iteration: ' num2str(kg)],'FontSize',24);
    try 
      [kg fm(45,35) fms(cur_ind(45,35) + cur_step(45,35))]
    end
    drawnow;
  end

end

