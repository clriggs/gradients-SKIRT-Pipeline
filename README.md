# NIHAO-SKIRT-Pipeline

Template respository for generating SKIRT radiative transfer mock observations from the NIHAO simulations.   

The methodology used here is described in  
["Panchromatic Simulated Galaxy Observations from the NIHAO Project"](https://iopscience.iop.org/article/10.3847/1538-4357/acf9f0 "The Astrophysical Journal")  

Analysis code for the resulting Fits files can be found at    
[https://github.com/ntf229/NIHAO-SKIRT-Catalog](https://github.com/ntf229/NIHAO-SKIRT-Catalog "NIHAO-SKIRT-Catalog") in analysis.py  

Starting from simulation outputs, code in the bin directory should be run in the following order:
1. makeParticles.py - create numpy arrays containing relevent simulation data
2. sampleOrientations.py - run RT for a sample of orientations for each galaxy uniformly on a sphere (r-band only)
3. calcAxisRatios.py - calculate axis ratios of the sampled orientations
4. selectOrientations.py - select orientations that have axis ratios which span each galaxy's range as evenly as possible
5. bestParamsRT.py - run full RT on selected orientations
6. storeBestAsFits.py - store mock photometry as Fits file
7. resolvedImagesFromFits.py - make images from Fits file with custom color scaling

