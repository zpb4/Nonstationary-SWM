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
### Executing program
The workflow below is configured to run from the file configuration when the Zenodo repository is unzipped and stored in a repository named 'data'
#### Hybrid SWM fitting
Numbering indicates order in which scripts must be run  
Runtimes (in parentheses at end) are estimated with parallelization where applicable on an HPC resource 

1) data_process.R: Processes raw hydrologic and state variable data for HYMOD 'process' model from .txt files in data repository (<1 min)
2) data_procss_sacsma.R: Processes raw hydrologic and state variable data for SAC-SMA 'truth' model from .txt files in data repository (<1 min)
3) hymod_error_pre-process: Data pre-processing routine for HYMOD data (<1 min)
4) sacsma_error_pre-process: Data pre-processing routine for SAC-SMA data (<1 min)
5) fit_model_hymod_cal-val.R: Fits hybrid SWM model to calibration-validation subsets of training data for HYMOD (<1 min)
6) lamc_synthetic-gen_hc.R: Generates specified number of synthetic ensemble samples saved in 'out' repository (6 hrs per 100 samples)
7) error_check-remove_hc.R: Synthetic ensemble post-processing script to remove/replace errant members (10 min per 100 samples)
8) syn-hefs_py-transfer_feather.R: Transfers R array ensemble output to .feather files for compatibility with Python (30 min per 100 samples)
9) ens_sample_npz-transfer.py: Arranges .feather files to Numpy array and zips individual synthetic ensembles for compatibility with EFO model (30 min per 100 samples)

#### Hybrid SWM generation

1) lamc_synthetic-gen_pre-hc.R: Generates specified number of pre-hindcast synthetic ensemble samples saved in 'out' repository (6 hrs per 100 samples)
2) error_check-remove_pre-hc.R: Synthetic ensemble post-processing script to remove/replace errant members, pre-hindcast period (10 min)
3) syn-hefs_py-transfer_feather_pre-hc.R: Transfers R array ensemble output to .feather files for compatibility with Python
4) ens_sample_npz-transfer.py: Arranges .feather files to Numpy array and zips individual synthetic ensembles for compatibility with EFO model

#### Functions
As above but for out-of-sample pre-hindcast generation

1) lamc_synthetic-gen_pre-hc.R: Generates specified number of pre-hindcast synthetic ensemble samples saved in 'out' repository (6 hrs per 100 samples)
2) error_check-remove_pre-hc.R: Synthetic ensemble post-processing script to remove/replace errant members, pre-hindcast period (10 min)
3) syn-hefs_py-transfer_feather_pre-hc.R: Transfers R array ensemble output to .feather files for compatibility with Python
4) ens_sample_npz-transfer.py: Arranges .feather files to Numpy array and zips individual synthetic ensembles for compatibility with EFO model

#### Plotting routines

- main_plot.R: Main plotting script for manuscript figures
- calc_ensemble_stats_top-10.R: Calculates cumulative ensemble statistics for top-10 inflow events
- calc_ensemble_stats_top-100.R: Calculates cumulative ensemble statistics for top-10 inflow events

#### Miscellaneous

#### Contact
Zach Brodeur, zpb4@cornell.edu
