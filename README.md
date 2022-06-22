# VirtualCT - 20220412 version
This project provides a virtual fan CT system with an elementary information for objects.
This is composed of two basic codes,

1. material based forward projection algorithm (MBFPA),   
2. Fltered back projection (FBP).

The former includes the simple beam-hardening correction algorithm as well as the signal noise model.

Plese cite the following papars for use.

Ref: Kai-Wen Li, et al., Physical density estimations of single- and dual-energy CT using material-based forward projection algorithm: a simulation study, Br J Radiol, 2021  
Ref: Kai-Wen Li, et al., kV-kV and kV-MV DECT based estimation of proton stopping power ratio - a simulation study, Physica Medica, v89, p182, 2021  
Ref: Daiyu Fujiwara, et al., Virtual computed-tomography system for deep-learning-based material decomposition, Physics in Medicine and Biology, 2022

**GPU version is not free, but can be distributed for some specific GPU boards. Contact us by e-mail. haga@tokushima-u.ac.jp


## 1: For use of material based forward projection algorithm (MBFPA) code
### 1-1: Preparation
Go to "Release_MBFPA_CPU_human/"
It is required to prepare an object including the elementary information and an X-ray spectrum before runing the code.
The elementary information of the ICRP110 human phantom and the digital Gammex phantom has been uploaded as the examples (Weight_input_ICRU110.raw and Weight_input_gammex.raw), both of which the datatype is the 16-bit unsigned, 512 $\times$ 512.
For human phantom, six major elements, H, C, N, O, P, and Ca, are included, whereas for Gammex phantom, eight major elements, H, C, N, O, P, Ca, Mg, and Si are included. We provided the creation code of digital Gammex phantom (phantom_creation_Gammex.cpp).
The 120 kV and 6 MV X-ray spectra were also provided within "Spectrum/" folder, which are simulated with a radiotherapy machine (ELEKTA Synergy). One can apply any energy spectrum with changing this file (first and second rows mean the first and last points in energy bin, and third row means the fraction for photons).

The information of the object and the spectrum must be indicated in the input txt file (IR_TOMO_input_virtual_projection_6MV.txt and IR_TOMO_input_virtual_projection_120kV.txt are the example).

The geometry can be varied by modifying parameters listed in physParamsTOMO.h, where
the original values simulate projections in the geometry of helical tomothrapy.

We note that the original code is for ICRP110 human phantom (6 elements). If one wants to simulate with Gammex phantom, "number of materials" and some lines of prior_weight_production.cpp need to be modified.

### 1-2: How to execute

0. Prepare elementary density image and incident X-ray spectrum  
      ex. Weight_input_ICRU110.raw, Spectrum/PhaseSpace_photon_120kV_0cm.txt  
1. make  
2. ./mvct.exe (input.txt)  
      ex. ./mvct.exe IR_TOMO_input_virtual_projection_120kV.txt        
3. output image is producted (as "reprojection_float.raw")  

The datatype of the projection image is a 32-bit real, 609 $\times$ 800 in the original geometory. 

### 1-3: To adjust a signal noise
The photon noise on the virtual detector is controlled by the parameter
"Photon noise in detector"
in the input txt file
The noise model is based on a Gaussian distribution. One can change it by modifing the function,
gaussiannoise0
in main.cpp.

### 1-4: To adjust a beam-hardening correction
We provided a beam-hardening correction model, which has been applied in Ref: Kai-Wen Li, et al., Physica Medica, v89, p182, 2021.
This is performed by running the code, "beam_hardening_correction.cpp".

Usage is here.

0. Prepare sinogram (projection image)  
      ex. reprojection_float.raw        
1. Compile as gcc beam_hardening_correction.cpp  
2. ./a.out "alpha" "beta"  
      ex. ./a.out 0.01 2.00 (for kV CT) or ./a.out 0.01 3.70 (for MV CT)        
3. output image is producted (as "reprojection_float_cor.raw")

Note that "alpha" and "beta" could depend on CT geometory(SDD, SID, Det.size etc.)/protocol(angle interval etc.) as well as the photon energy and phantom size.
Therefore, it is recommended that these parameters are optimized by a homogeneity check with a water phantom or something like that.


## 2: For use of Fltered back projection (FBP) code
### 2-1: Preparation
Go to "Release_FBP_CPU/".
It is required for CT reconstruction to prepare the projection data (sinogram).
An example is "reprojection_float.raw" or "reprojection_float_cor.raw".
The information of the sinogram must be indicated in the input txt file as,
"IR_TOMO_input_virtual_projection_6MV.txt" for MVCT reconstrcution.
Note that the same geometry as that used in the MBFPA should be employed there.

### 2-2: How to execute

0. Prepare sinogram (projection image)  
      ex. reprojection_float.raw from "Release_MBFPA_CPU_human/"  
1. make  
2. ./mvct.exe (input.txt)  
      ex. ./mvct.exe TOMO_input_virtual_projection_human.txt  
3. output image is producted (as "FBP_virtual_projection_512x512_human.raw")  

The datatype of the reconstructed image is a 32-bit real, 512 $\times$ 512. 
