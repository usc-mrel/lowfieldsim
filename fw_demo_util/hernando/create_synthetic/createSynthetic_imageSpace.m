%% Function: createSynthetic_imageSpace
%%
%% Description: create synthetic chemical shift-encoded dataset, with arbitrary species and echo times. 
%%
%% Some features:
%%    - Image-space
%%    - Single-R2*
%%    - Accepts multi-peak species
%%    - Accepts multiple species
%%    - Accepts a field map in Hz
%%
%% Input arguments:
%%   - imDataParams0.TEs: echo times (in seconds)
%%   - imDataParams0.fieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
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
%%   - trueParams.species(ii).amps: true water/fat images, size [nx,ny,nz,ncoils] 
%%   - trueParams.r2starmap: R2* map (in s^{-1}, size [nx,ny,nz])
%%   - trueParams.fieldmap: field map (in Hz, size [nx,ny,nz])
%%
%%
%% Output:
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
%%
%% Author: Diego Hernando
%% Date created: August 19, 2011
%% Date last modified: November 10, 2011
%%


function imDataParams = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams )

gyro = 42.58;
imDataParams = imDataParams0;
t = imDataParams0.TE;
[sx,sy,sz,C] = size(trueParams.species(1).amps);
N = length(t);


try
  fieldStrength = imDataParams0.FieldStrength;
catch
  fieldStrength = 1.5;
end


  
try 
  r2starmap = trueParams.r2starmap;
catch 
  r2starmap = zeros(sx,sy);
end
  
try 
  fieldmap = trueParams.fieldmap;
catch 
  fieldmap = zeros(sx,sy);
end
  
imDataParams.images = zeros(sx,sy,sz,C,N);
for ks=1:length(trueParams.species)
  amps = trueParams.species(ks).amps;

  try 
    freqs = gyro*fieldStrength*algoParams.species(ks).frequency;
  catch 
    freqs = 0;
  end
   
  try
    relAmps = algoParams.species(ks).relAmps;
  catch
    relAmps = ones(size(freqs))/length(freqs);
  end
  
  s = zeros(1,1,1,1,N);
  fieldAndR2starEffect = zeros(sx,sy,sz,C,N);
  for kt=1:N
    s(1,1,1,1,kt) = sum(relAmps(:).*exp(j*2*pi*freqs(:)*t(kt)));
    fieldAndR2starEffect(:,:,:,:,kt) = exp(-t(kt)*repmat(r2starmap,[1 1 1 C]) + j*2*pi*t(kt)*repmat(fieldmap,[1 1 1 C]));
  end
  
  imDataParams.images = imDataParams.images + repmat(amps,[1 1 1 1 N]).*repmat(s,[sx sy sz C 1]).*fieldAndR2starEffect;
end


%%   - imDataParams.PrecessionIsClockwise (1 = fat has positive frequency; -1 = fat has negative frequency)
try
  imDataParams.PrecessionIsClockwise = imDataParams0.PrecessionIsClockwise;
  if imDataParams.PrecessionIsClockwise <= 0 
    imDataParams.PrecessionIsClockwise == -1;
  end
catch 
  imDataParams.PrecessionIsClockwise = -1;
end

% If precession is clockwise, make everything rotate in the opposite direction
if imDataParams.PrecessionIsClockwise < 0
  imDataParams.images = conj(imDataParams.images);
end

