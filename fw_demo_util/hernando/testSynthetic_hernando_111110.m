%% testSynthetic_hernando_111110
%%
%% Test fat-water algorithms on synthetic data
%%
%% Author: Diego Hernando
%% Date created: August 19, 2011
%% Date last modified: February 29, 2012

% Add to matlab path
BASEPATH = './';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);

%% Create some synthetic "true" data, and set acquisition parameters
sx = 192;sy=sx;
imtest = phantom(sx);
threshold = 0.8;
N = 6; TEinit = 1.2e-3; dTE = 2e-3;
TE = TEinit + [0:N-1]*dTE;
trueParams.species(1).amps = imtest.*(imtest<threshold); % Water
trueParams.species(2).amps = imtest.*(imtest>=threshold); % Fat

[X,Y] = meshgrid(linspace(-1,1,sx),linspace(-1,1,sy));
trueParams.fieldmap = 20*randn*ones(sx,sy) + 40*randn*X + 40*randn*Y + 100*randn*X.^2 + 100*randn*Y.^2  + 100*randn*X.*Y.^2 + 100*randn*X.^3  + 100*randn*Y.^3;
trueParams.r2starmap = 80 - 40*imtest;

imDataParams0.TE = TE;
imDataParams0.FieldStrength = 1.5;
imDataParams0.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];


% Simulate data
imDataParams = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams );

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 120]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut 
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );

% $$$ %% Recon -- mixed fit for phase error correction
% $$$ % Initialize mixed fitting to graph cut solution
% $$$ algoParams.fieldmap = outParams.fieldmap;
% $$$ algoParams.r2starmap = outParams.r2starmap;
% $$$ algoParams.NUM_MAGN = 1;
% $$$ algoParams.THRESHOLD = 0.04;
% $$$ algoParams.range_r2star = [0 200];
% $$$ 
% $$$ % Do mixed fitting
% $$$ outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams )



