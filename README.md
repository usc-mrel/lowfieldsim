<div align="center">
<a href="http://www.opensourceimaging.org/project/minimum-field-strength-simulator/" target="_blank">
<img src="https://usc-mrel.github.io/logo/lfs_logo_small.jpg" width="300"></img>
</a>
</div>

MATLAB Scripts for Low SNR Simulation at Low Field Strengths
============================================================
Under certain assumptions, low-field MRI acquisitions can be simulated from 
high-field data and minimum field strength requirements can be determined. This 
package presents a simple framework for simulating low-field MRI acquisitions, 
which can be used to predict the minimum B0 field strength necessary for MRI 
techniques. This framework may be particularly useful for the evaluation of 
de-noising and constrained reconstruction techniques, and the possibilities of 
translating them to less expensive low-field scanners. 

(c) *Weiyi Chen*, *Ziyue Wu*, *Krishna Nayak*, May 2016.

[Magnetic Resonance Engineering Laboratory](https://mrel.usc.edu)

**University of Southern California**

Code Structure
--------------
### Main Function 
**lowfieldgen.m**

	function [ k_low ] = lowfieldgen( inParam )
	% LOWFIELDGEN simulates low field noise
	% See details inside the m-file

### Demo
We recommend start using this package by running these 2 demos:

1. Upper Airway Gridding Reconstruction:	**demo_airway.m**  
   This example shows gridding reconstruction on simulated low field data, based on data 
   acquired using Golden angle radial FLASH at 3T.
	
2. Fat-water Separation:	**demo_fatwater.m**  
   This example shows fat-water separation on simulated low field data, based on data 
   acquired using product IDEAL GRE sequence at 3T.
		
Assumptions
-----------
1. **Body noise dominance** : 
	We assume that body thermal noise is the dominant noise source at all field 
	strengths under investigation (0.1 – 3.0 T).
2. **Consistent B<sub>1</sub><sup>+</sup> field** :
	We assume that the uniformity of RF transmission is consistent across field 
	strengths.
3. **Consistent B<sub>1</sub><sup>-</sup> field** :
	We assume that the receiver coils have the same geometry and noise covariance 
	at different field strengths.
4. **Consistent B<sub>0</sub> homogeneity** :
	We assume the same off-resonance in parts-per-million (ppm) at different field
	 strengths. This results in less off-resonance in Hertz at lower field.
5. **Single species dominance** :
	We use a single global relaxation correction function to account for the 
	signal change at different field strengths.
6. **Steady state acquisition** :
	If the signals are not acquired at steady state, the magnetization relaxation 
	will be determined not only by the sequence parameters but also by the initial
	 state.

Publications
--------------
1. Wu Z, Chen W, Nayak KS (2016) Minimum Field Strength Simulator for Proton Density 
   Weighted MRI. PLoS ONE 11(5): e0154711. doi:10.1371/journal.pone.0154711. 
   [JNRL](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0154711)
2. Wu Z, Chen W, Nayak KS (2016) Low-Field Simulation and Minimum Field Strength 
   Requirements. ISMRM Workshop on Data Sampling and Image Reconstruction, 
   Sedona, Arizona.
   [TALK](http://www.ismrm.org/workshops/Data16/program.htm) 

Acknowledgment
--------------
Fat-water separation is implemented using the graph cut field-map estimation 
method from the 
[ISMRM fat-water toolbox](http://ismrm.org/workshops/FatWater12/data.htm)

License
--------------
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0"
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/"
property="dct:title">Minimum Field Strength Simulator</span> by <span xmlns:cc="http://creativecommons.org/ns#"
property="cc:attributionName">Weiyi Chen, Ziyue Wu, Krishna Nayak</span> is licensed under a <a rel="license"
href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
