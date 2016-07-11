% Low Field Simulation demo - upper airway
% Data acquired using Golden angle radial FLASH at 3T
% Gridding reconstruction

% Ziyue Wu, Aug 2014
% Weiyi Chen, May 2016
% University of Southern California.

%% House keeping
clear; close all;

load upperairway3T_short.mat

%% Parameters for the simulations
inParam.B_high   = 3;       % acquired data field stength
inParam.B_low    = 0.5;     % target field strength
inParam.tissue   = 'muscle';
inParam.sequence = 'GradientEcho';
inParam.theta    = 5;
inParam.BW_high  = 62.5;
inParam.BW_low   = 62.5;
inParam.TR_high  = 4.6;
inParam.TR_low   = 4.6;
inParam.n_cov    = N_cov;
inParam.TE_high = 2.6;
inParam.TE_low = 2.6;

%% Add noise to simulate low field strength
inParam.k_high = reshape(k_rad, [cv.Nkx cv.Nky 1 1 cv.Ncoil]);

% simulate low field data
k_low = lowfieldgen(inParam);
k_rad_low = reshape(k_low, [cv.Nkx cv.Nky cv.Ncoil]);

% keep the high field accquired data just for comparison
k_rad_high = reshape(inParam.k_high, [cv.Nkx cv.Nky cv.Ncoil]);

%% Recon parameters:
windowsize = 21; %temporal window size
% # of views b/t 2 temporal frames, equal to windowsize means no
% viewsharing
viewinterval = 21;
maxFrame = floor(size(k_rad_low,2)/windowsize) -1;
fprintf('Maximum frames to be recontructed is %d.\n', maxFrame)
init_frame = 11; % choose from 0 ~ [max_frame-1]
if init_frame < 0 || init_frame > maxFrame
    error('wrong init_frame!');
end

final_frame =11;
if final_frame < init_frame || final_frame > maxFrame
    error('wrong final_frame!');
end
% if set 0:  do recon from 0 to (max_frame-1)
tot_frame = final_frame - init_frame + 1;
fprintf('Reconstructing from frame #%d to #%d.\n', init_frame, final_frame)

frameshow = 1; % # of frames to be displayed, must <= tot_frame
if frameshow > tot_frame
    error('wrong frameshow!');
end
os = 2; % oversampling factor in gridding
kernalsize = 2*os; % *os to be in oversampled grid

%% Gridding
tot_spoke = size(k_rad_low, 2);
if tot_frame > 0
    % do nothing, use specified tot_frame
else
    if windowsize >= viewinterval
        tot_frame = floor((tot_spoke-windowsize)/viewinterval);
    else
        tot_frame = floor(tot_spoke/(viewinterval));
    end
end

img_low = zeros(cv.Nkx*os,cv.Nkx*os,cv.Ncoil,tot_frame);
img_high = zeros(cv.Nkx*os,cv.Nkx*os,cv.Ncoil,tot_frame);
imgsize = cv.Nkx;

for tempframe = init_frame: final_frame
    cur_T = T(:,(tempframe*viewinterval + 1) ...
        : (tempframe*viewinterval + windowsize));
    cur_W = W(:,(tempframe*viewinterval + 1) ...
        : (tempframe*viewinterval + windowsize));
    cur_k_rad_low = k_rad_low(:,(tempframe*viewinterval + 1) ...
        : (tempframe*viewinterval + windowsize),:);
    cur_k_rad_high = k_rad_high(:,(tempframe*viewinterval + 1) ...
        : (tempframe*viewinterval + windowsize),:);
    
    for coil = 1:cv.Ncoil
        img_low(:,:,coil,tempframe - init_frame +1) ...
            = gridkb(cur_k_rad_low(:,:,coil),cur_T,cur_W,...
            imgsize,os,kernalsize,'image');
        img_high(:,:,coil,tempframe - init_frame +1) ...
            = gridkb(cur_k_rad_high(:,:,coil),cur_T,cur_W,...
            imgsize,os,kernalsize,'image');
    end
    fprintf('Reconstructing frame #%d.\n', tempframe);
end

%% display images
fprintf('Displaying from frame #%d to #%d.\n', ...
    init_frame, init_frame + frameshow -1)
% low field
img_low_combined = sqrt(sum(abs(img_low.^2), 3));
img_low_combined = permute(img_low_combined, [ 1 2 4 3]);
img_low_combined = img_low_combined(129:(128+256) , 129:(128+256),:);

% accquired
img_high_combined = sqrt(sum(abs(img_high.^2), 3));
img_high_combined = permute(img_high_combined, [ 1 2 4 3]);
img_high_combined = img_high_combined(129:(128+256) , 129:(128+256),:);

if(frameshow > 0)
    for tempframe = 1: frameshow
        figure; 
        
        % accquired
        subplot(121),
        imshow(abs(img_high_combined(:,:,tempframe)),[]);
        set(gca,'FontSize',18);
        title('Accquired @ 3T');
        
        % low field
        subplot(122),
        imshow(abs(img_low_combined(:,:,tempframe)),[]);
        set(gca,'FontSize',18);
        title(['Simulated @ ',num2str(inParam.B_low),'T']);
    end
end