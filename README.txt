Deconvolution pipeline
======================
Figure S5 presents the details of the pipeline.

-First step: align_PSF_DataStack.m

Empirical 3D system PSF is aligned (along the z-axis) with the raw data stack. 
This is achieved by performing deconvolution (Richardson-Lucy) of a small number of 
z-slices (typically separated by 100 µm) of data with a set of 2D PSFs sampled at 
different depths (typically separated by 10 µm). The resulting deconvolved images 
are analyzed (manually or automatically) for sharpness to determine the global 
z-axis alignment of the system PSF and the raw data stack.

-Second step: deconv_Time_Series_Data.m
All the z-slices of the image stack are deconvolved using 2D PSFs sampled from system 
PSF at correspondingly aligned z-positions. For time lapse datasets, PSF and stack 
alignment is calculated using the first time point data, which is then used to 
deconvolve all the time points



Segmentation
============
nuclear_Segmentation.m
script used for cellular segmentation



Delta F over F
==============
gen_Baseline_For_DFOF.m : generates baseline F as an average over entire recording duration
gen_DFOF.m : generate delta F over F


PCA and ICA analyses for identifying synchrony (Figure 7)
=================================================
Synchrony_ICA_PCA_analysisCode.m : All the code used to generate Figure 7 is
calculate_iterative_noise.m : function to calculate the noise level cutoff for dF/F traces


Utility scripts
===============
proj_In_Time.m : generate maximum intensity, average intensity and standard deviation over a 
specified recording duration
gen_XZ_YZ_MIProj.m : generate XZ and YZ maximum intensity projections

empirical_psf_analysis.py
Python scripts used for FWHM calculations of the empirical PSF  (Figure 1C-top plot)
