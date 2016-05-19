% Function: decomposeGivenFieldMapAndDampings
%
% Description: estimate water/fat images given the nonlinear parameters
% 
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%
%  - fieldmap: the estimated B0 field map
%  - r2starWater: the estimated water R2* map
%  - r2starFat: the estimated fat R2* map
%
% Returns: 
%  - amps: the amplitudes for all chemical species and coils
%  - remerror: fit error norm
%
% Author: Diego Hernando
% Date created: 
% Date last modified: August 18, 2011

function [amps,remerror] = decomposeGivenFieldMapAndDampings( imDataParams,algoParams,fieldmap,r2starWater,r2starFat )

gyro = 42.58;



try
  precessionIsClockwise = imDataParams.PrecessionIsClockwise;
catch
  precessionIsClockwise = 1;
end

try 
  ampW = algoParams.species(1).relAmps;
catch
  ampW = 1.0
end

  
% If precession is clockwise (positive fat frequency) simply conjugate data
if precessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;

sx = size(images,1);
sy = size(images,2);
N = size(images,5);
C = size(images,4);

relAmps = reshape(relAmps,1,[]);


B1 = zeros(N,2);
B = zeros(N,2);
for n=1:N
  B1(n,:) = [ampW*exp(j*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(j*2*pi*deltaF(2:end)*t(n)))];
end

remerror = zeros(sx,sy);
for kx =1:sx
  for ky=1:sy
    s = reshape( squeeze(images(kx,ky,:,:,:)), [C N]).';

    B(:,1) = B1(:,1).*exp(j*2*pi*fieldmap(kx,ky)*t(:) - r2starWater(kx,ky)*t(:));
    B(:,2) = B1(:,2).*exp(j*2*pi*fieldmap(kx,ky)*t(:) - r2starFat(kx,ky)*t(:));

    amps(kx,ky,:,:) = B\s;

    if nargout > 1
      remerror(kx,ky) = norm(s - B*squeeze(amps(kx,ky,:,:)),'fro');
    end
  end
end


