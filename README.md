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
1. Upper Airway Gridding Reconstruction:	**demo_airway.m**
2. Fat-water Separation:	**demo_fatwater.m**
		
Assumptions
-----------
1. **Body noise dominance** : 
	We assume that body thermal noise is the dominant noise source at all field 
	strengths under investigation (0.1 – 3.0 T).
2. **Consistent B1^+ field** :
	We assume that the uniformity of RF transmission is consistent across field 
	strengths.
3. **Consistent B1^- field** :
	We assume that the receiver coils have the same geometry and noise covariance 
	at different field strengths.
4. **Consistent B0 homogeneity** :
	We assume the same off-resonance in parts-per-million (ppm) at different field
	 strengths. This results in less off-resonance in Hertz at lower field.
5. **Single species dominance** :
	We use a single global relaxation correction function to account for the 
	signal change at different field strengths.
6. **Steady state acquisition** :
	If the signals are not acquired at steady state, the magnetization relaxation 
	will be determined not only by the sequence parameters but also by the initial
	 state.

Acknowledgment
--------------
Fat-water separation is implemented using the graph cut field-map estimation 
method from the 
[ISMRM fat-water toolbox](http://ismrm.org/workshops/FatWater12/data.htm)