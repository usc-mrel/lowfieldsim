% Function name: createExpansionGraphVARPRO_fast
%
% Description: Create the adjacency graph for the fieldmap update
% problem. Basically, solving a minimum-cut problem on this graph is
% equivalent to finding the best among the 2^Q neighbors of the
% current fieldmap estimate.
%
% Input arguments:
%   - residual: the 3D residual array (of size NUM_FMS x SX x SY)
%   - dfm: (scalar) the discretization step for the fieldmap: dfm = fms(n) - fms(n-1)
%   - lambda: (SX x SY) the map of regularization parameters
%   - size_clique: the size of the MRF clique (typically 1)
%   - cur_ind: (map of SX x SY positive integers) the current field map (as a map of indices)
%   - step: (map of SX x SY integers, all of the same sign) the current jump move
%
% Output arguments:
%   - A: (sparse matrix of size (SX*SY+2)x(SX*SY+2) nonnegative integers) the desired graph adjacency matrix
%   - Aunsc: the desired graph adjacency matrix before scaling/rounding 
%
%Author: Diego Hernando
%Date created: March 18, 2008
%Date last update: October 21, 2009
%--------------------------------------------------------------------------

function [A,Aunsc] = createExpansionGraphVARPRO_fast( residual, dfm, lambda, size_clique, cur_ind, step )

% Allocate space for the adjacency matrix A (which completely describes the graph)
sx = size(residual,2);
sy = size(residual,3);
L = size(residual,1);
s = sx*sy;
num_nodes = s + 2;
num_edges = num_nodes*(2 + (size_clique+1)^2);

sA = [num_nodes, num_nodes];

% Auxiliary arrays
offset = [0:(s-1)]'*L; 
step_ind = cur_ind + step;

% Scaling
if strcmp(computer,'GLNX86')
  maxA = 1e5;
else
  maxA = 1e6;
end

valsh = zeros(1,s);
valsv = zeros(s,1);

% Insert smoothness-related edges in graph
factor = lambda*dfm^2;
x = 1:sx;
y = 1:sy;
[Y,X] = meshgrid(y,x);

allIndCross = [];
allValsCross = [];
for dx = -size_clique:size_clique
  for dy = -size_clique:size_clique
    dist = sqrt(dx^2+dy^2);
    if dist>0
      validmapi = X+dx>=1 & X+dx<=sx & Y+dy>=1 & Y+dy<=sy;
      validmapj = X-dx>=1 & X-dx<=sx & Y-dy>=1 & Y-dy<=sy;
      
      
      curfactor = min(factor(validmapi),factor(validmapj)); 
      
      a = curfactor.*(1./dist*(cur_ind(validmapi)-cur_ind(validmapj)).^2);
      b = curfactor.*(1./dist*(cur_ind(validmapi)-step_ind(validmapj)).^2);
      c = curfactor.*(1./dist*(step_ind(validmapi)-cur_ind(validmapj)).^2);
      d = curfactor.*(1./dist*(step_ind(validmapi)-step_ind(validmapj)).^2);

      
      temp = zeros(1,s);
      temp(validmapi) = max(0,c-a);
      valsh = valsh + temp;

      temp = zeros(s,1); 
      temp(validmapi) = max(0,a-c);
      valsv = valsv + temp;
      
      temp = zeros(1,s); 
      temp(validmapj) = max(0,d-c);
      valsh = valsh + temp;
      
      temp = zeros(s,1); 
      temp(validmapj) = max(0,c-d);
      valsv = valsv + temp;
    
      S = 1:s;
      Sh = 1 + (S(validmapi(:)));
      Sv = 1 + (S(validmapj(:)));
      indcross = sub2ind(sA, Sh,Sv );
      temp = b+c-a-d;

      allIndCross = [allIndCross,indcross];
      allValsCross = [allValsCross,temp(:)'];
    end
  end
end



ind1 = sub2ind(sA, ones(1,num_nodes-2),2:(num_nodes-1) );

ind2 = sub2ind(sA,2:(num_nodes-1) , num_nodes + zeros(1,num_nodes-2));

% Insert residual-related edges in graph 
temp0 = residual(cur_ind(:)+offset(:));
valid_ind = (step_ind(:)>=1 & step_ind(:)<=L);
temp1 = zeros(s,1);
temp1(valid_ind~=0) = residual(step_ind(valid_ind~=0)+offset(valid_ind~=0));
curmaxA = max(max(temp0),max([valsh.';valsv;allValsCross.']));
infty = curmaxA;
temp1(valid_ind==0) = infty;

indAll = [ind1 ind2 allIndCross];
valuesAll = [valsh + reshape(max(temp1-temp0,0),1,s), valsv.' + reshape(max(0,temp0-temp1),1,s), allValsCross];

[indAllSort,sortIndex] = sort(indAll);

% This way of creating my matrix seems faster
valsort = valuesAll(sortIndex);
[xind,yind] = ind2sub([num_nodes num_nodes], indAllSort);
A = sparse(xind,yind,valsort,num_nodes,num_nodes,length(valsort));

% Rescale and round adjacency matrix
if nargout>1
  Aunsc = A; 
end

A = round(A*maxA/curmaxA);
A(A<0) = 0; % Remove residual entries < 0
