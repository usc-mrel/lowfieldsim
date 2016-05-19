function [m, p] = gridkb(d,k,w,n,osf,wg,opt)

% function m = gridkb(d,k,w,n,osf,kw,opt)
%
%     d -- k-space data
%     k -- k-trajectory, scaled -0.5 to 0.5
%     w -- k-space weighting
%     n -- image size (m will be osf*n X osf*n)
%     osf -- oversampling factor (usually between 1 and 2)
%     wg -- full kernel width in oversampled grid samples (usually 3 to 7)
%     opt -- 'k-space', 'image', defaults it 'image' if not specified
%
%     m -- gridded k-space data
%     p -- gridding kernel, optional
%
%  Uses optimum Kaiser-Bessel window for a given
%    oversampling factor and kernel size
%  Now uses Phil's numbers
%


%  Written by John Pauly, 2003, 2005, 2007, 2011
%  (c)Board of Trustees, Leland Stanford Junior University
%  Modified by Ziyue Wu September 2012. --- add a constant to reduce deapodization at the edge of FOV.

if nargin < 7,
    opt = 'image';
end
    
% convert to single column
d = d(:);
k = k(:);
w = w(:);

% width of the kernel on the original grid
kw = wg/osf;

% preweight
dw = d.*w;

% compute kernel, assume e1 is 0.001, assuming nearest neighbor
kosf = floor(0.91/(osf*1e-3));

% half width in oversampled grid units
kwidth = osf*kw/2;

% beta from the Beatty paper
beta = pi*sqrt((kw*(osf-0.5)).^2-0.8);

% compute kernel
om = [0:kosf*kwidth]/(kosf*kwidth);
p = besseli(0,beta*sqrt(1-om.*om));
p = p./p(1);
% last sample is zero so we can use min() below for samples bigger than kwidth
p(end) = 0;

% convert k-space samples to matrix indices
nx = (n*osf/2+1) + osf*n*real(k);
ny = (n*osf/2+1) + osf*n*imag(k);

m = zeros(osf*n,osf*n);

% loop over samples in kernel at grid spacing
for lx = -kwidth:kwidth,
  for ly = -kwidth:kwidth,

    % find nearest samples
    nxt = round(nx+lx);
    nyt = round(ny+ly);

    % seperable kernel value
    kkx = min(round(kosf*abs(nx-nxt)+1), floor(kosf*kwidth)+1);
    kwx = p(kkx);
    kky = min(round(kosf*abs(ny-nyt)+1), floor(kosf*kwidth)+1);
    kwy = p(kky);

    % if data falls outside matrix, put it at the edge, zero out below
    nxt = max(nxt,1); nxt = min(nxt,osf*n);
    nyt = max(nyt,1); nyt = min(nyt,osf*n);

    % accumulate gridded data
    m = m+sparse(nxt,nyt,dw.*kwx'.*kwy',osf*n,osf*n);
  end;
end;

% zero out data at edges, which is probably due to data outside mtx
m(:,1) = 0; m(:,osf*n) = 0;
m(1,:) = 0; m(osf*n,:) = 0;

% stop here, if we just want the k-space data
if strcmp(opt,'k-space') return; end;

im = ifftshift(ifft2(ifftshift(m)));


% compute deappodization function
x = [-osf*n/2:osf*n/2-1]/(n);
sqa = sqrt(pi*pi*kw*kw*x.*x-beta*beta);
dax = sin(sqa)./(sqa);
% normalize by DC value
dax = dax/dax(osf*n/2);
% make it a 2D array
da = dax'*dax;

% deappodize
% im = im./da;
im = im./(da + 1); % add a constant to reduce deapodization at the edge of FOV, ZW.

%return the result
m = im;
end
