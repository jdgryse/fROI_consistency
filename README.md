# Welcome

This repository accompanies the paper "Improving consistency of functional regions of interest in fMRI"

1)  This repository holds both code and data files that can be used across platforms (linux, Windows, macOS) to perform all methods mentioned in the original paper. The results of the paper are obtained with the following operating system: macOS Mojave, version 10.14.4.
2)  - All .R files contain code that should be run with the statistal program R. In order to run this code, the user should install R or RStudio. The version of R used for the results in the paper is R 3.5.3. Once R is installed, the code can readily be run in the console, with the exception for the paths where data are stored locally, which should be adjusted accordingly,  and additional R packages that may not already be installed locally.
    - All .sh files are shell scripts (Unix language), with code that can be run in a Bash shell. In order for most
       commands to work, FSL should be installed on your computer. The results in this paper were obtained with FSL 5.0.9.

# single_subject: Analyses with data from Gonzalez-Castillo (2012)

## figures

See README.txt in the folder.

## toy_data

A toy data set in order to perform the consistency analysis performed in the original study, as well as an overall brain mask for this “simulated individual”, since the data described in Section 2.2 of the original paper are not freely available. These are all NIFTI files that can be read into R or FSL or the fMRI data analysis program of your choice. 

## scripts

### analysis.R

In this file, NHST, the ABT and mLR method are performed on every run under different parameter configurations, which is described in Section 2.2 of the paper. The result of this file are thresholded SPMs. The resulting SPMs were then used to calculate the different consistency measures mentioned in Section 2.2.

### consistency_coherence.R

Cohen's kappa for NHST, ABT and the mLR is computed based on the combined thresholded maps created in create_coherencemaps.R

### consistency_maitra.R

The Maitra similarity index is computed on every possible pair of two runs.

### consistency_voxels.R

The number of detected voxels in each of the runs is computed as well as the standard deviation of the number of detected voxels across all runs.

### create_coherenemaps.R

Create the combined thresholded maps of all runs in order to compute Cohen's kappa in consistency.coherence.R

### design.

These files contain information concerning the design of the fixed effects analysis in FSL, needed for the flameo command in real_es_groundtruth.sh and do not need to be adjusted in order to work for the toy data set.

### make_conditions.R

This file makes a table with all parameter configurations needed to reproduce the parameter combinations for the ABT and mLR method mentioned in Section 2.2 in the paper.  The resulting tables are conditions_abt.txt and conditions_mLR.txt. These tables are used in consistency_coherence.R, consistency_maitra.R, consistency_voxels.R, create_coherencemaps.R and TDR.R

### real_designCrossVal.R

Here, the design. files mentioned above are made.

### real_es_groundtruth.sh

Here, the ground truth of effect sizes for each step of the cross-validation is constructed. The ground truth of effect sizes is then used in analysis.R.

### real_writefiles.R

In this file, the overall brain mask for the “simulated individual” is constructed. This is then used in  analysis.R.

### TDR.R

The truth detection rate, described in Section 2.2 of the paper, is computed in order to evaluate sensitivity of NHST, the ABT and the mLR method.


# HCP: Analyses with data from the Human Connectome Project


## copes

The estimated contrast estimates of both runs from the first level analysis for 10 of the 50 genetically unrelated subjects analyzed in the paper. This is described in Section 2.3.

## varcopes

The estimated (standard error)^2 of the contrast estimates of both runs from the first level analysis for 10 of the 50 genetically unrelated subjects analyzed in the paper. This is described in Section 2.3.

## hybrid_mask

This folder contains all masks and code to obtain these masks described in Section 2.3 of the paper.

### *_fusiform.nii.gz

The anatomical masks that were obtained using the Harvard-Oxford probabilistic atlas in FSL and were combined in a bilateral mask of the fusiform gyrus.

### *_hemisphere.nii.gz

The anatomical masks for the left and right hemisphere. These were used to get the right FFA.

### FFA.nii.gz

Anatomical bilateral mask of the fusiform gyrus created using the *_fusiform.nii.gz masks and hybrid_mask.R

### FFA_functional.nii

Functional mask of the FFA based the forward-inference meta-analysis conducted on Neurosynth.

### FFA_hybrid_right.nii.gz

Mask representing the overlap between FFA.nii.gz, right_hemisphere.nii.gz and FFA_functional.nii. This was created in hybrid_mask.R and used in all analyses in the folder ../scripts.

### hybrid_mask.R

File to combine the separate *_fusiform.nii.gz masks into the bilateral anatomical fusiform gyrus mask FFA.nii.gz. The mask representing the overlap between FFA.nii.gz, right_hemisphere.nii.gz and FFA_functional.nii is created as well.

## scripts

This folder contains all code to perform the procedure described in Section 2.3 and to obtain the Figures in Section 3.2.1 and 3.2.2

### analysis_hcp.R

The contrast estimates and their standard errors from ../copes and ../varcopes are used to peform NHST, ABT and mLR for different parameter configurations on both runs per subject. This is done for both the general and summary rFFA fROI. Consistency measures calculated here are the Maitra similarity index for both the general and summary rFFA fROI as well as euclidean distance between peaks for the general rFFA fROI between the two runs for each subject. This is described in Section 2.3 of the paper. The resulting arrays are then used to create the Figures in Sections 3.2.1 and 3.2.2. Functions that are used here can be found in functions_hcp.R

### figures_hcp_cond.R

Using the output of analysis_hcp.R, the Figures of Sections 3.2.1 and 3.2.2 are created here.

### functions_hcp.R

Functions that are needed for analysis_hcp.R: function for calculating the Maitra similarity index, functions for calculating the euclidean distance between peaks defined with NHST, ABT and the mLR method and functions to define the summary rFFA fROI based on the peak defined with NHST, ABT and the mLR method.

### make_conditions.R

This file makes a table with all parameter configurations needed to reproduce the parameter combinations for the ABT and mLR method mentioned in Section 2.3 in the paper.  The resulting tables are conditions_abt_hcp.txt and conditions_mLR_hcp.txt. These tables are used in analysis.R




