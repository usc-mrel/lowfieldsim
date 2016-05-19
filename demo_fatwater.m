% Low Field Simulation demo - fat-water separation
% Data acquired using IDEAL GRE at 3T
%
% Fat-water separation based on Magn Reson Med. 2010 Jan;63(1):79-90
% see also fw_i2cm1i_3pluspoint_hernando_graphcut, test_hernando_111110
%
% (c) Aug. 2014, Weiyi(Wayne) Chen, University of Southern California

%% House keeping
clear; clc; close all;
load 'fat-water@3T-3echo.mat';
addpath( genpath( './fw_demo_util' ) );

B0_low = 0.3;
%% Low SNR Simulation

inParam.B_high   = 3;
inParam.B_low    = B0_low;
inParam.tissue   = 'liver';
inParam.sequence = 'GradientEcho';
inParam.theta    = 3;
inParam.BW_high  = 62.5;
inParam.BW_low   = inParam.BW_high * B0_low/3;
inParam.TR_high  = 9;
inParam.TR_low   = TE(end) * 1000 * 3/B0_low + 1/inParam.BW_low/2 + 5.608;
inParam.n_cov    = n_cov;

% simulation for each TE
k_low = zeros(size(k_high));
for t = 1:3
    inParam.k_high   = k_high(:,:,:,t,:);
    inParam.TE_high  = TE(t) * 1000;                % TE in ms
    % Lengthen TE to reach at the same level of off-resonance
    inParam.TE_low   = TE(t) * 1000 * 3/B0_low;     
    k_low(:,:,:,t,:) = lowfieldgen(inParam);
end

% fat-water toolbox requires data in format [Nkx Nky Nkz Ncoil NTE]
k_low = permute(k_low, [1 2 3 5 4]);

%% Assamble imDataParams structure
recon = sqrt( size(k_low,1) * size(k_low,2) )...
        * ifftshift(ifft2(fftshift(k_low)));
imDataParams.images = recon;
imDataParams.FieldStrength = 3;
imDataParams.TE = TE;
imDataParams.PrecessionIsClockwise = 1;

%% Separation parameters setup
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency =[-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood 
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation 
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent 
algoParams.LMAP_POWER = 2; % Spatially-varying regularization
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Field inhomogeneities using a graph cut algorithm. 
% Magn Reson Med. 2010 Jan;63(1):79-90.
tic
  outParams ...
      = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc
fat_frac = rot90( computeFF(outParams), 3 );
water = rot90( outParams.species(1).amps, 3 );
fat = rot90( outParams.species(2).amps, 3 );

%% Output
figure;
subplot(131),imshow(abs(water),[]); title('Water'); 
subplot(132),imshow(abs(fat),[]); title('Fat'); 
subplot(133),imshow(fat_frac,[0 100]);
title('Fat Fraction');

