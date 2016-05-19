Author: Diego Hernando
Date: February 14, 2012

This folder contains functions for image-space fat-water separation
and fat quantification. The main functions contained are:

* graphcut/fw_i2cm1i_3pluspoint_hernando_graphcut.m
Fat-water separation using a regularized formulation and discretized
graphcut algorithm

(Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)

* descent/fw_i2cm0i_3plusploint_hernando_optimtransfer.m
Fat-water separation using a regularized formulation and continuous
descent "optimization-transfer" algorithm. Can be run automatically 
at the end of the graphcut function to remove field map quantization
effects. Based on:

(Huh W, Fessler JA, Samsonov AA. Water-fat decomposition with regularized field map. 
In: Proceedings of the 16th Annual Meeting of ISMRM, Toronto, Canada, 2008. p. 1382.)

* mixed_fitting/fw_i2xm1c_3pluspoint_hernando_mixedfit.m
Fat-water separation voxel by voxel, including correction for phase
errors by "mixed fitting"

(Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging 
using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)

* create_synthetic/createSynthetic_imageSpace.m
Creation of synthetic chemical shift-encoded datasets, based on a
ground truth (fat, water, R2*, field map), model (fat peaks), and
acquisition parameters (echo times, field strength)

* common/coilCombine.m
Coil combination for image sequences, based on:

(Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of phased
array MR imagery. Magn Reson Med 2000;43:682-690)

* common/coilCombine3D.m
Coil combination for 3D image sequences, simply calling the 2D function described above, for each slice. 




The top-level folder contains two test scripts showing how to call the above functions: 

* test_hernando_111110.m (test graphcut + descent and mixed fitting on
  randomly-chosen dataset from Peter Kellman)

* testSynthetic_hernando_111110.m (test synthetic dataset creation,
  based on Shepp-Logan phantom, and graphcut + descent separation)



Some issues: 

* Compatibility of matlab_bgl package for solving the graphcut
  formulation is not guaranteed with all platforms (see below)

* So far, has been tested succesfully on the following platforms:

-- Linux, 64bit (R2009b) 

-- Linux, 32bit (R2009b) 

-- Windows, 64bit (R2009a, max_flow may not be working properly)

-- MacOS, 64 bit
   