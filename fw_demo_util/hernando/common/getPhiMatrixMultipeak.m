% Function: getPhiMatrixMultipeak
% 
% Returns the Phi matrix of complex exponentials from the VARPRO
% formulation of Dixon imaging
% 
% Arguments:
%   - deltafs: chemical shifts of the different chemical species
%   - relAmps: relative amplitudes of different fat peaks
%   - t: echo times used
% 
% Returns: 
% 
%   - Phi: the Phi matrix, Phi_{i,j} =  exp(j 2 \pi deltaf_j t_i )
%
%
% Author: Diego Hernando
% Date created: Aug 27, 2008
% Last modified: Aug 27, 2008
%
function Phi = getPhiMatrixMultipeak( deltafs,relAmps, t )

[DF,T] = meshgrid( deltafs,t );
[A,T2] = meshgrid( relAmps,t );


Phi1 = exp(j*2*pi*T.*DF);
Phi = [Phi1(:,1) , sum(Phi1(:,2:end).*A,2)];
