%% Function: checkParamsAndSetDefaults_optimtransfer
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

function [validParams,algoParams2] = checkParamsAndSetDefaults_optimtransfer( imDataParams,algoParams )

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


%% Initial guess for R2* map
try 
  algoParams2.r2starmap = algoParams.r2starmap;
catch
  algoParams2.r2starmap = zeros(size(imDataParams.images(:,:,1,1,1,1)));
end

%% Initial guess for field map
try 
  algoParams2.fieldmap = algoParams.fieldmap;
catch
  algoParams2.fieldmap = zeros(size(imDataParams.images(:,:,1,1,1,1)));
end

%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
try
  algoParams2.size_clique = algoParams.size_clique;
catch 
  algoParams2.size_clique = 1;
end

%%   - algoParams.OT_ITERS = 30; % Number of optimization transfer iterations
try
  algoParams2.OT_ITERS = algoParams.OT_ITERS;
catch 
  algoParams2.OT_ITERS = 30;
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


