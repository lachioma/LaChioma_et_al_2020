# Code accompanying La Chioma et al., 2020

This repository contains code used to perform data analyses and to reproduce figures of the article  [Disparity sensitivity and binocular integration in mouse visual cortex areas. La Chioma et al. 2020](https://www.biorxiv.org/content/10.1101/2020.05.07.083329v1).

The main script `analysis_code_LaChioma2020.m` is organized in separate sections that allow analysis and reprodution of individual figures. 
Each section starts with the data files that need to be loaded for reproducing that specific figure.


## Where to get the data
Associated data files are located at https://edmond.mpdl.mpg.de/imeji/collection/_chMNO_mc1EQ5Gu5

Place the data `.mat` files in the folder \Vars, for loading directly from the main script (`analysis_code_LaChioma2020.m`).


## How data files are organized

The data files contain processed data after extraction of fluorescence time courses from two-photon images.

N.B. Not all data files are needed for each analysis. Each section of the main script `analysis_code_LaChioma2020.m` starts with the data files that are needed for reproducing that specific figure.

- Files var_DGD_V1.mat, var_DGD_LM.mat, var_DGD_RL.mat contain data obtained using dichoptic gratings (DGD stands for Drifting Gratings Dichoptic), for  experiments  in each area (V1, LM, RL).
\
Each var_DGD file contains a structure variable named aDGD, with the following fields:
  * ROIs:
    * ver: raw calcium traces for each ROI
    * dFoF: dF/F for each ROI
    * nTrials: actual nr. of trials for each stimulus
    * reps: nr. of trials set for the experiment (if exp is interrupted, nTrials will be lower than reps)
    * ExpID: experiment ID for each ROI
    * There are many other fields with computed data (mean, max, std, etc.). The field name should be self-explanatory.

  * Param: some parameters used for processing fluorescence data, for each experiment.

  * StimSettings: settings of dichoptic gratings, for each experiment.

  * Info:  info about the experiments, for each experiment.

  * FP: values of two-photon images (pixel size, um, etc.)

  * ROI: values of the ROI segmentation, for each ROI


- Files aD_RDS_V1.mat, aD_RDS_LM.mat, aD_RDS_RL.mat contain data obtained using random dot stimuli (RDS), for experiments in each area (V1, LM, RL). 
The RDS variables are structured in a similar way to DGD variables described above.

- Files aD_RPS_V1.mat, aD_RPS_LM.mat, aD_RPS_RL.mat contain data with mapping of receptive field elevation, for experiments in each area (V1, LM, RL). 
The RPS variables are structured in a similar way to DGD variables described above.

- File Disp_DGD contains disparity selectivity index (DI) values and related quantifications for each experiment and at different levels of stringency of response inclusion criteria.

- File DIvar_DGD contains disparity selectivity index (DI) values.

- Files 'Fit' contain tuning curves and fitted curves. Fits are computed only for responsive ROIs (see Methods; thr_std=4, thr_rep_fraction=0.5).

- File NC_DGD contains noise correlation analysis.

- File NN_DGD contains spatial organization analysis.

- File OD_DG contains quantifications of ocular dominance.

- File R_DGD contains quantifications of dF/F, facilitation and suppression indexes.

- Folder SVM discrimination contains population decoding analysis with support vector machine (SVM).

- File TC_A__DGD contains tuning curves obtained with DGD and distributions of preferred disparieties.

N.B. Throughout codes and data, the abbreviations RDC and RDS are used interchangeably to refer to random dot stimuli. 
Note that the random dot stimuli used in La Chioma et al. 2020 (here referred to as random dot correlograms RDC) are the same as in [La Chioma et al. 2019](http://dx.doi.org/10.1016/j.cub.2019.07.037) (random dot stereograms RDS).
The correct name for these stimuli is random dot correlograms (thanks to the reviewer for pointing this out!).
