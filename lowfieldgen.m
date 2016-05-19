function [ k_low ] = lowfieldgen( inParam )
%LOWFIELDGEN_TEST simulates low field noise
%|  [ k_low ] = lowfieldgen_test( inparam )
%|
%|  Output:
%|      k_low:      simulated kspace data at low field in format 
%|                  [Nkx Nky Nkz Nt Ncoil]
%|
%|  Input:
%|      inparam:    input parameter structure, details below:
%|
%|      .k_high:        kspace data aquired at high field in format 
%|                      [Nkx Nky Nkz Nt Ncoil]
%|      .B_high(T):     B0 field strength at which data was acquired
%|                      choose from (0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 3) 
%|                      or other values (need to specify .T1_high)
%|      .B_low(T):      simulated B0 field strength
%|                      choose from (0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 3) 
%|                      or other values (need to specify .T1_low)
%|      .tissue:        'muscle'
%|                      'kidney'
%|                      'white matter'
%|                      'gray matter'
%|                      'liver'
%|                      'fat'
%|                      'other': need to specify .T1_low .T1_high .T2
%|      .sequence:      currently support:
%|                      'SpinEcho'
%|                      'GradientEcho'
%|                      'bSSFP'
%|                      'InversionRecovery'
%|      .TR_high(ms):   TR of kspace data
%|      .TR_low(ms):    TR of simulated low field data
%|      .TE_high(ms):   TE of kspace data
%|      .TE_low(ms):    TE of simulated low field data
%|      .BW_high(kHz):  Readout bandwidth of kspace data
%|      .BW_low(kHz):   Readout bandwidth of simulated low field data
%|      .theta(degree): flip angle
%|      .n_cov:         noise covariance matrix [Ncoil Ncoil],
%|                      if a single number is entered, a diagonal matrix 
%|                      will be used
%|      (optional):
%|      .T1_high(s):    required only if T1 values if tissue type is
%|                      'other' or B_high is NOT from 
%|                      (0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 3)
%|      .T1_low(s):     required only if T1 values if tissue type is
%|                      'other' or B_low is NOT from 
%|                      (0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 3)
%|      .T2(s):         required only if T2 values if tissue type is
%|                      'other'
%
% (c) written by Ziyue Wu, Feburary 2014.
% (c) modified by Weiyi Chen, August 2014.
% University of Southern California
% https://mrel.usc.edu

%% T1 & T2 Correction Table     REF: Principles of MRI. D.G. Nishimura
B0_set = [0.1 0.2 0.3 0.4 0.5 1 1.5 3];

tissue_type = {'muscle'; 'kidney'; 'white matter'; 'gray matter';...
    'liver'; 'fat'; 'other'};

T1_muscle =  [0.28 0.37 0.45 0.50 0.55 0.73 0.87 1.42];
T1_kidney = [0.34 0.39 0.44 0.47 0.50 0.58 0.65 1.19];
T1_wm = [0.31 0.39 0.45 0.49 0.53 0.68 0.78 1.08];
T1_gm = [0.40 0.49 0.55 0.61 0.65 0.82 0.92 1.82];
T1_liver = [0.18 0.23 0.27 0.30 0.32 0.43 0.49 0.81];
T1_fat = [0.17 0.18 0.20 0.21 0.22 0.24 0.26 0.30];

T2_muscle = 0.047;
T2_kidney = 0.058;
T2_wm = 0.092;
T2_gm = 0.1;
T2_liver = 0.043;
T2_fat = 0.085;

%% Passing parameters

% Covariance matrix
[Nkx, Nky, Nkz, Nt, Ncoil] = size(inParam.k_high);
if size(inParam.n_cov) == 1
    n_cov = inParam.n_cov * eye(Ncoil);
else
    n_cov = inParam.n_cov;
end

% T1 at high field
if ~ismember('T1_high',fieldnames(inParam))
    if ~any(B0_set == inParam.B_high) ||...
            ~ismember(inParam.tissue, tissue_type)
        error ( 'specify T1 value at high field strength');
    else
        [~,ind] = ismember(inParam.B_high,B0_set);
        switch (inParam.tissue)
            case 'muscle'
                T1_high = T1_muscle(ind);
            case 'kidney'
                T1_high = T1_kidney(ind);
            case 'white matter'
                T1_high = T1_wm(ind);
            case 'gray matter'
                T1_high = T1_gm(ind);
            case 'liver'
                T1_high = T1_liver(ind);
            case 'fat'
                T1_high = T1_fat(ind);
            case 'other'
                error ( 'specify T1 value at high field strength');
        end
    end
else
    T1_high = inParam.T1_high;
end

% T1 at low field
if ~ismember('T1_low',fieldnames(inParam))
    if ~any(B0_set == inParam.B_low) ||...
            ~ismember(inParam.tissue, tissue_type)
        error ( 'specify T1 value at low field strength');
    else
        [~,ind] = ismember(inParam.B_low,B0_set);
        switch (inParam.tissue)
            case 'muscle'
                T1_low = T1_muscle(ind);
            case 'kidney'
                T1_low = T1_kidney(ind);
            case 'white matter'
                T1_low = T1_wm(ind);
            case 'gray matter'
                T1_low = T1_gm(ind);
            case 'liver'
                T1_low = T1_liver(ind);
            case 'fat'
                T1_low = T1_fat(ind);
            case 'other'
                error ( 'specify T1 value at low field strength');
        end
    end
else
    T1_low = inParam.T1_low;
end

% T2
if ~ismember('T2',fieldnames(inParam))
    if ~ismember(inParam.tissue, tissue_type)
        error ( 'specify T2 value');
    else
        switch (inParam.tissue)
            case 'muscle'
                T2 = T2_muscle;
            case 'kidney'
                T2 = T2_kidney;
            case 'white matter'
                T2 = T2_wm;
            case 'gray matter'
                T2 = T2_gm;
            case 'liver'
                T2 = T2_liver;
            case 'fat'
                T2 = T2_fat;
            case 'other'
                error ( 'specify T2 value');
        end
    end
else
    T2 = inParam.T2;
end

% TE
TE_high = inParam.TE_high / 1000;       %ms --> sec
TE_low = inParam.TE_low / 1000;         %ms --> sec

% TR
TR_high = inParam.TR_high / 1000;       %ms --> sec
TR_low = inParam.TR_low / 1000;         %ms --> sec

% Flip angle
theta = inParam.theta * pi / 180;       %degree --> rad

% Readout bandwidth
BW_high = inParam.BW_high;
BW_low = inParam.BW_low;

%% Signal Scaling

E1_h = exp( -TR_high / T1_high );
E1_l = exp( -TR_low / T1_low );
E2_h = exp( -TE_high / T2 );
E2_l = exp( -TE_low / T2 );

a = inParam.B_low / inParam.B_high;

switch inParam.sequence
    case 'SpinEcho'       
        fx = ( ( 1 - E1_l ) / ( 1 - E1_l * cos(theta) ) ) ...
            / ( ( 1 - E1_h ) / ( 1 - E1_h * cos(theta) ) ) ...
            * ( E2_l / E2_h );
        scaleS = fx * (a^2);
    case 'GradientEcho'
        fx = ( ( 1 - E1_l ) / ( 1 - E1_l * cos(theta) ) ) ...
            / ( ( 1 - E1_h ) / ( 1 - E1_h * cos(theta) ) ) ...
            * ( E2_l / E2_h );
        scaleS = fx * (a^2);
    case 'bSSFP'
        
    case 'InversionRecovery'
        
end

%% Noise Scaling

b = BW_low / BW_high;

scaleN = sqrt(a^2*b - (a^4)*(fx^2));
[v,d] = eig(n_cov);
% shuffle random seeds so that randn doesn't give the same # each time!
rng('shuffle'); 
noise = scaleN * (v*sqrt(d))*randn(Ncoil,Nkx*Nky*Nkz*Nt);
noise = permute(noise, [2 1]);
noise = reshape(noise, [Nkx Nky Nkz Nt Ncoil]);

%% Output

k_low = inParam.k_high * scaleS + noise;



end