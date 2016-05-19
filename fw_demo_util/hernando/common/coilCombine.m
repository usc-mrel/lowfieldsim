% Function: coilCombine
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,nz,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nz,1,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function im2 = coilCombine( im1 )

% Let's make the coil dimension the fourth one and the TE the third
im1 = permute(im1,[1 2 5 4 3]);

% Get image dimensions and set filter size
[sx,sy,N,C] = size(im1);
filtsize = 7;

% Initialize
im2 = zeros(sx,sy,1,1,N);
Rs = zeros(sx,sy,C,C);

% Get correlation matrices
for kc1=1:C
  for kc2=1:C
    for kn=1:N
      Rs(:,:,kc1,kc2) = Rs(:,:,kc1,kc2) + filter2(ones(filtsize),im1(:,:,kn,kc1).*conj(im1(:,:,kn,kc2)),'same');
    end
  end
end

% Compute and apply filter at each voxel
for kx=1:sx
  for ky=1:sy
% $$$     [U,S] = eig(squeeze(Rs(kx,ky,:,:)));
% $$$     s = diag(S);
% $$$     [maxval,maxind] = max(abs(s));
% $$$     myfilt = U(:,maxind);    
% $$$     im2(kx,ky,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);

    % Change suggested by Mark Bydder
    [U,S] = svd(squeeze(Rs(kx,ky,:,:)));
    myfilt = U(:,1); 
    im2(kx,ky,1,1,:) = myfilt'*reshape(squeeze(im1(kx,ky,:,:)).',[C N]);
  end
end

% In case the input data are single
if strcmp(class(im1),'single')
  im2 = single(im2);
end
