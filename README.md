# Nonstationary-SWM
Code repository to support WRR manuscript 'A hybrid, non-stationary Stochastic Watershed Model (SWM) for uncertain hydrologic projections under climate change'   
Submitted 4 April 2023

## Description
The code below supports the data processing, model fitting, generation, and plotting routines to support the aforementioned manuscript.
## Getting started
### Dependencies
Releases of this code are stored in the following Zenodo repository: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7702354.svg)](https://doi.org/10.5281/zenodo.7702354)
### Installing
Requires following R packages:
* ranger
* fGarch
* doParallel
* optimx
### Executing program
The workflow below is configured to run from scripts in the 'src' folder to fit and then simulate from the hybrid SWM. Additional scripts for benchmarking, helper functions, and figure pre-processing are included in the 'Miscellaneous Functions' block and plotting scripts for primary manuscript (ms) and supporting information (si) figures are included in the 'Plotting' block.
#### Hybrid SWM fitting
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) src/data_process.R: Processes raw hydrologic and state variable data for HYMOD 'process' model from .txt files in data repository (<1 min)
2) fit_model_hymod.R: Fits hybrid SWM model to calibration-validation subsets of training data for HYMOD (5 min)
6) fit_model_hymod_cal-val_benchmark.R: Fits static SWM model to calibration-validation subsets of training data for HYMOD (<1 min)
7) fit_model_sacsma_cal-val_skip-samp.R: Fits hybrid SWM model to calibration-validation subsets of training data for SAC-SMA for 'skip-sample' approach outlined in manuscript (5 min)
7) fit_model_sacsma_cal-val_skip-samp.R: Fits hybrid SWM model to calibration-validation subsets of training data for SAC-SMA for 'split-sample' approach outlined in manuscript (5 min)

#### Hybrid SWM simulation

1) src/syn_gen_hymod.R: Generates specified number of hybrid SWM samples for HYMOD in Test scenario; saved in 'out' repository (6 hr per 1000 samples)
2) src/syn_gen_hymod_4c.R: Generates specified number of hybrid SWM samples for HYMOD in Test+4C scenario; saved in 'out' repository (6 hr per 1000 samples)
3) syn_gen_sacsma.R: Generates specified number of hybrid SWM samples for SAC-SMA; saved in 'out' repository (6 hr per 1000 samples)
4) syn_gen_sacsma_split-samp.R: Generates specified number of hybrid SWM samples for SAC-SMA in with 'split-sample' approach; saved in 'out' repository (6 hr per 1000 samples)

#### Functions

1) GL_maineqs_mv.R: Required functions to support fitting and generation from dynamic residual model
2) GL_subeqs.R: Subequations for 'GL_maineqs_mv.R'

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu
