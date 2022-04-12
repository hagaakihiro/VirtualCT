# VirtualCT - 20220412 version
This project provides a virtual fan CT system with an elementary information for objects.
This is composed of two basic codes,
1. material based forward projection algorithm (MBFPA),
2. Fltered back projection (FBP).
The former includes the simple beam-hardening correction algorithm as well as the signal noise model.


## 1: For use of material based forward projection algorithm (MBFPA) code
### 1-1: Preparation
It is required to prepare an object including the elementary information and an X-ray spectrum before runing the code.
The elementary information of the ICRP110 human phantom and the digital Gammex phantom has been uploaded as the examples (Weight_input_ICRU110.raw and Weight_input_gammex.raw), both of which the datatype is the 16-bit unsigned, 512 \times 512.
For human phantom, six major elements, H, C, N, O, P, and Ca, are included, whereas for Gammex phantom, eight major elements, H, C, N, O, P, Ca, Mg, and Si are included. We provided the creation code of digital Gammex phantom (phantom_creation_Gammex.cpp).
The 120 kV and 6 MV X-ray spectra were also provided within "Spectrum/" folder, which are simulated with a radiotherapy machine (ELEKTA Synergy). One can apply any energy spectrum with changing this file (first and second rows mean the first and last points in energy bin, and third row means the fraction for photons).

The information of the object and the spectrum must be indicated in the input txt file (IR_TOMO_input_virtual_projection_6MV.txt and IR_TOMO_input_virtual_projection_120kV.txt are the example).

We note that the original code is for ICRP110 human phantom (6 elements). If one wants to simulate with Gammex phantom, "number of elements" and some lines of prior_weight_production.cpp need to be modified. We are sorry for this inconvenient.

### 1-2: How to execute

### 1-3: To adjust a signal noise

### 1-4: To adjust a beam-hardening correction


## 2: For use of Fltered back projection (FBP) code
### 2-1: Preparation

### 2-2: How to execute
