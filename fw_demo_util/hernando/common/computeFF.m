function ff = computeFF( outParams, ffParams )

if nargin<2
  ffParams.noise_bias_correction = 1;
end

        
curf = outParams.species(2).amps;
curw = outParams.species(1).amps;


denom = (abs(curf) + abs(curw));
denom2 = denom;
denom2(denom==0) = 1; % To avoid divide-by-zero issues
ff = 100*abs(curf)./denom2;


if ffParams.noise_bias_correction>0
  fatregions = ff>50;
  watregions = ff<=50;
  denom2 = abs(curf + curw);
  denom2(denom==0) = 1; % To avoid divide-by-zero issues
  ff(watregions) = 100 - 100*abs(curw(watregions))./denom2(watregions);
  ff(fatregions) = 100*abs(curf(fatregions))./denom2(fatregions);
end