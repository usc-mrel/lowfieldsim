% Function name: getQuadraticApprox
%
% Description: Estimate the second derivative of the residual (as a
% function of estimated field map), at the minimizer at each voxel
%
% Input arguments:
%   - residual: the 3D residual array (of size NUM_FMS x SX x SY)
%   - dfm: (scalar) the discretization step for the fieldmap: dfm = fms(n) - fms(n-1)
%
% Output arguments:
%   - d2: the second derivatives of the residual (at the minimizer at each voxel)
%
%Author: Diego Hernando
%Date created: June 9, 2009
%Date last update: June 9, 2009
%--------------------------------------------------------------------------

function d2 = getQuadraticApprox( residual, dfm )

[NUM_FMS,sx,sy] = size(residual);
resoffset = [0:(sx*sy-1)]'*NUM_FMS;

[minres,iminres] = min(residual(10:end-9,:,:),[],1);
iminres = squeeze(iminres + 9);

d2 = (residual((iminres(:)+1) + resoffset) + residual((iminres(:)-1) + resoffset) - 2*residual((iminres(:)) + resoffset)  )/dfm^2;

d2 = reshape(d2,[sx sy]);

