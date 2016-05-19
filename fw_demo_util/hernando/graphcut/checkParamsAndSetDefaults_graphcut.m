%% Function: checkParamsAndSetDefaults_graphcut
%%
%% Description: Check validity of input parameters and set defaults for unspecified parameters
%%
%% Input:
%%   - imDataParams: TEs, images and field strength
%%   - algoParams: algorithm parameters
%%
%% Output:
%%   - validParams: binary variable (0 if parameters are not valid for this algorithm)
%%   - algoParams2: "completed" algorithm parameter structure (after inserting defaults for unspecified parameters)
%% 
%%
%% Author: Diego Hernando
%% Date created: August 19, 2011
%% Date last modified: November 10, 2011
%%

function [validParams,algoParams2] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams )

imDataParams2 = imDataParams;
algoParams2 = algoParams;
validParams = 1;

% Start by checking validity of provided data and recon parameters
if size(imDataParams,3) > 1
  disp('ERROR: 2D recon -- please format input data as array of size SX x SY x 1 X nCoils X nTE')
  validParams = 0;
end

if length(algoParams.species) > 2
  disp('ERROR: Water=fat recon -- use a multi-species function to separate more than 2 chemical species')
  validParams = 0;
end

if length(imDataParams.TE) < 3
  disp('ERROR: 3+ point recon -- please use a different recon for acquisitions with fewer than 3 TEs')
  validParams = 0;
end

%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
try
  algoParams2.size_clique = algoParams.size_clique;
catch 
  algoParams2.size_clique = 1;
end

%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
try
  algoParams2.range_r2star = algoParams.range_r2star;
catch 
  algoParams2.range_r2star = [0 0];
end

%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
try
  algoParams2.NUM_R2STARS = algoParams.NUM_R2STARS;
catch 
  algoParams2.NUM_R2STARS = 1;
end

%%   - algoParams.range_fm = [-400 400]; % Range of field map values
try
  algoParams2.range_fm = algoParams.range_fm;
catch 
  algoParams2.range_fm = [-400 400]; 
end

%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
try
  algoParams2.NUM_FMS = algoParams.NUM_FMS;
catch 
  algoParams2.NUM_FMS = 301;
end

%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
try
  algoParams2.NUM_ITERS = algoParams.NUM_ITERS;
catch 
  algoParams2.NUM_ITERS = 40;
end

%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
try
  algoParams2.SUBSAMPLE = algoParams.SUBSAMPLE;
catch 
  algoParams2.SUBSAMPLE = 1;
end

%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
try
  algoParams2.DO_OT = algoParams.DO_OT;
catch 
  algoParams2.DO_OT = 0;
end

%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
try
  algoParams2.LMAP_POWER = algoParams.LMAP_POWER;
catch 
  algoParams2.LMAP_POWER = 2; 
end

%%   - algoParams.lambda = 0.05; % Regularization parameter
try
  algoParams2.lambda = algoParams.lambda;
catch 
  algoParams2.lambda = 0.05;
end

%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
try
  algoParams2.LMAP_EXTRA = algoParams.LMAP_EXTRA;
catch 
  algoParams2.LMAP_EXTRA = zeros(size(imDataParams.images(:,:,1,1,1)));
end

%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
try
  algoParams2.TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
catch 
  algoParams2.TRY_PERIODIC_RESIDUAL = 0;
end

%%   - imDataParams.PrecessionIsClockwise (1 = fat has positive frequency; -1 = fat has negative frequency)
try
  imDataParams2.PrecessionIsClockwise = imDataParams.PrecessionIsClockwise;
  if imDataParams2.PrecessionIsClockwise <= 0 
    imDataParams2.PrecessionIsClockwise == -1;
  end
catch 
  imDataParams2.PrecessionIsClockwise = -1;
end


