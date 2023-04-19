# Nonstationary-SWM
Code repository to support WRR manuscript 'A hybrid, non-stationary Stochastic Watershed Model (SWM) for uncertain hydrologic projections under climate change'   
Submitted 4 April 2023

## Description
The code below supports the data processing, model fitting, generation, and plotting routines to support the aforementioned manuscript.
## Getting started
### Dependencies
Raw data to support this code in in the following Zenodo repository: https://doi.org/10.5281/zenodo.7689054
### Installing
Requires following R packages:
* ranger
* fGarch
*doParallel
### Executing program
The workflow below is configured to run from the file configuration when the Zenodo repository is unzipped and stored in a repository named 'data'
#### Hybrid SWM fitting
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) data_process.R: Processes raw hydrologic and state variable data for HYMOD 'process' model from .txt files in data repository (<1 min)
2) data_procss_sacsma.R: Processes raw hydrologic and state variable data for SAC-SMA 'truth' model from .txt files in data repository (<1 min)
3) hymod_error_pre-process: Data pre-processing routine for HYMOD data (<1 min)
4) sacsma_error_pre-process: Data pre-processing routine for SAC-SMA data (<1 min)
5) fit_model_hymod_cal-val.R: Fits hybrid SWM model to calibration-validation subsets of training data for HYMOD (5 min)
6) fit_model_hymod_cal-val_benchmark.R: Fits static SWM model to calibration-validation subsets of training data for HYMOD (<1 min)
7) fit_model_sacsma_cal-val_skip-samp.R: Fits hybrid SWM model to calibration-validation subsets of training data for SAC-SMA for 'skip-sample' approach outlined in manuscript (5 min)
7) fit_model_sacsma_cal-val_skip-samp.R: Fits hybrid SWM model to calibration-validation subsets of training data for SAC-SMA for 'split-sample' approach outlined in manuscript (5 min)

#### Hybrid SWM generation

1) syn_gen_hymod_tst.R: Generates specified number of hybrid SWM samples for HYMOD in Test scenario; saved in 'out' repository (1 hr per 1000 samples)
2) syn_gen_hymod_4c.R: Generates specified number of hybrid SWM samples for HYMOD in Test+4C scenario; saved in 'out' repository (1 hr per 1000 samples)
3) syn_gen_sacsma_skip-samp.R: Generates specified number of hybrid SWM samples for SAC-SMA in with 'skip-sample' approach; saved in 'out' repository (1 hr per 1000 samples)
4) syn_gen_sacsma_split-samp.R: Generates specified number of hybrid SWM samples for SAC-SMA in with 'split-sample' approach; saved in 'out' repository (1 hr per 1000 samples)

#### Functions

1) GL_maineqs_mv.R: Required functions to support fitting and generation from dynamic residual model
2) GL_subeqs.R: Subequations for 'GL_maineqs_mv.R'

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures
- calc_ensemble_stats_top-10.R: Calculates cumulative ensemble statistics for top-10 inflow events
- calc_ensemble_stats_top-100.R: Calculates cumulative ensemble statistics for top-10 inflow events

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu
