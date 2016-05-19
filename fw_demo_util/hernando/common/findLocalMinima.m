% Function name: findLocalMinima
%
% Description: Find the local minima of the VARPRO residual at each
% voxel, and return them in a 3D array.
%
% Input arguments:
%   - residual: the 3D residual array (of size NUM_FMS x SX x SY)
%   - threshold: a threshold for the signal at each voxel to be even considered 
%     (voxels with signal lower than the threshold will not be assigned any local 
%     minima, as these will be meaningless anyway)
%   - masksignal: alternatively, one can provide a binary map specifying which 
%     voxels have non-negligible signal
%
% Output arguments:
%   - masksignal: the masksignal (same as input, or the one computed inside 
%     function if not provided)
%   - resLocalMinima: the main output of this function: a 3D array giving, at 
%     each voxel, the indices where the local minima of the residual are located.
%   - numMinimaPerVoxel: SXxSY map specifying how many local minima were found 
%     in the residual at each voxel
%
%Author: Diego Hernando
%Date created: March 18, 2008
%Date last update: March 18, 2008
%--------------------------------------------------------------------------

function [masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima( residual, threshold, masksignal )


% Get the dimensions
L = size(residual,1);
sx = size(residual,2);
sy = size(residual,3);

% Take finite differences in residual
dres = diff(residual,1,1);

maxres = squeeze(max(residual,[],1));
minres = squeeze(min(residual,[],1));


% Find voxels with signal
if nargin < 3
  sumres = sqrt(squeeze(sum(residual,1)));
  sumres = sumres/max(max(sumres));
  masksignal = sumres>threshold;
end

% Loop through all voxels (this is MATLAB-inefficient and should be done in vectorized way)
resLocalMinima = zeros(1,sx,sy);
numMinimaPerVoxel = zeros(sx,sy,1);
for kx=1:sx
  for ky=1:sy
    if masksignal(kx,ky) > 0
      minres = min(residual(:,kx,ky));
      maxres = max(residual(:,kx,ky));

      temp = [0;squeeze(dres(:,kx,ky))];
      temp = temp<0 & circshift(temp,-1)>0 & residual(:,kx,ky)<minres+0.3*(maxres-minres);
      
      resLocalMinima(1:sum(temp),kx,ky) = find(temp);
      numMinimaPerVoxel(kx,ky) = sum(temp);
    end
  end
end

% $$$ x = 1:sx;
% $$$ y = 1:sy;
% $$$ [Y,X] = meshgrid(y,x);
% $$$ 
% $$$ for kx=1:sx
% $$$   for ky=1:sy
% $$$     if masksignal(kx,ky) == 0
% $$$ 
% $$$       noise_vals = 5:10:size(residual,1);
% $$$ 
% $$$       resLocalMinima(1:length(noise_vals),kx,ky) = noise_vals';
% $$$       numMinimaPerVoxel(kx,ky) = length(noise_vals);
% $$$     end
% $$$   end
% $$$ end
% $$$ 
