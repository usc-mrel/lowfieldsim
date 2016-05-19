% Function: coilCombine3D
%
% Description: combine multi-coil image sequences
%
% Based on: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of
% phased array MR imagery. Magn Reson Med 2000;43:682-690
% 
% Parameters:
% im1: the multi-coil images (size [nx,ny,ncoils,nTE])
%
% Returns: 
% im2: the coil-combined images (size [nx,ny,nTE])
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: February 7, 2012

function im2 = coilCombine3D( im1 )


[sx,sy,sz,C,N] = size(im1);

% Maintain the data type (e.g., single, double) of the input data
ims = zeros([sx,sy,sz,1,N],class(im1));
for kz=1:sz
  im2(:,:,kz,1,:) = coilCombine(im1(:,:,kz,:,:));
end
