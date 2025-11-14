# Schrödinger-Brain
We present a framework for augmenting brain dynamics into a complex-valued field using Hamilton's equations. We then employ a data-driven strategy to model brain dynamics using Schrödinger-like equation. The entire data-driven modeling (iincluding linear and nonlinear) and analysis of complex-valued brain network dynamics is grounded in multi-modal data, including resting-state and task fMRI, diffusion MRI, and auxiliary benchmarks.
## Overview
This repository contains the custom commercial, open-source, and in-house code used to analyse the data in our study. Most analyses were performed in **MATLAB R2018b**, with additional models and statistics implemented in **Python 3.9.15** and **R 4.2.3**.

The code is organised to allow:

- preprocessing of neuroimaging time series at the voxel or region-of-interest (ROI) level,
- estimation of coupling parameters and other network measures using schrodinger-based model, 
- statistical modelling of age-related and other covariate effects,
- generation of figures and tables reported in the manuscript.


## System Requirements
### Hardware requirements
- At least 16 GB RAM for region-level analyses; 32 GB or more is recommended for voxel-level analyses.
- Multi-core CPU is beneficial for parallel processing.
- Several tens of GB of free disk space depending on which datasets you use.

### Software requirements

#### OS requirements
This package is supported for Windows.
- Windows 64-bit

#### Matlab Dependencies
- **MATLAB R2018b** or later

Most core analyses and custom functions are implemented in MATLAB.

#### Python Dependencies
- **Python 3.9.15** (3.9+ recommended)
- Key packages and versions used in the study:

```text
numpy       1.24.0
torch       2.7.1
scipy       1.8.1
matplotlib  3.6.2
```
#### R Dependencies
- **R 4.2.3**
- Packages:
```text
mgcv  1.8-41   # generalised additive models for nonlinear age–connectivity associations
```
## Installation Guide
### Clone the repository 
```bash
git clone https://github.com/ShirleyZhang111/Schrodinger-Brain.git
cd Schrodinger-Brain
```

### Python environment
Create a dedicated environment and install the required packages, for example using `conda`:

```bash
conda create -n schrodinger-brain python=3.9
conda activate schrodinger-brain

pip install numpy==1.24.0 scipy==1.8.1 matplotlib==3.6.2 torch
# or, if a requirements file is provided:
# pip install -r requirements.txt
```
### R packages
In R:

```r
install.packages("mgcv")
```

### MATLAB path 
Within MATLAB, add the repository (and subfolders containing functions) to your path:

```matlab
addpath(genpath(pwd));
savepath;
```
## Data
The repository does **not** distribute any raw neuroimaging or third‑party datasets. To reproduce the results, you will need to obtain the data from the original sources and respect their usage policies.
### Human Connectome Project (HCP)
- Young Adults (HCP-YA) S1200 release  
- Lifespan Development and Aging (HCP Lifespan)
Access is available through the Human Connectome Project platforms; please follow the official instructions on the HCP website.

### UK Biobank imaging data
Neuroimaging data from UK Biobank are available by application via the UK Biobank access management system (project application ID as stated in the manuscript).

### Traffic benchmark data
Traffic-related time series used as a non-neural benchmark are downloaded from the **UCI Machine Learning Repository**.

## Quick Demo 
The following minimal examples demonstrate how to run a basic analysis after preparing your data. Here, we provide data for an individual subject from the HCP, HCPex, and UKB datasets, respectively, and execute the model. Note that the actual script names may vary depending on the final repository structure; please adjust paths and filenames accordingly.

### Python: Hamilotin Model

```python
# Example script: run demo_Hamilton.py
```

### MATLAB: Linear/Nonlinear Schrödienger-Llike model
```matlab
% Example script1 : run demo_Linear_Shcrodinger.m
% data-driven modeling using nonlinear Schröienger-like model

% Load rs-fMRI signals (available datasets: HCP-voxel)
load('TC_HCP_voxel.mat');
% Set model configuration parameters
field = 'complex';
inCfg = struct('field',field, 'T', 101:400, 'TC',TC); % Configuration structure: use time points 101-400 to compute model parameters
% Calculate parameters for linear Schrodinger-like model
[Q,G,lambda] = calc_coefficients_linear(inCfg);
```
```matlab
% Example script2 : run demo_Linear_Shcrodinger.m
% data-driven modeling using nonlinear Schröienger-like model

% Load rs-fMRI signals (available datasets: HCP, HCPex, HCP-voxel, UKB)
load('TC_HCP_379.mat');
% Set configuration parameters for nonlinear model
cfg.TC = TC;
cfg.is_largescale = 0;
cfg.maxit = 2000;
cfg.tol = 1e-12;
cfg.field = 'complex';
cfg.T = 101:1100;
cfg.mu = 0.1*length(cfg.T);
% Estimate Model Parameters from rs-fMRI data using nonlinear model
[coeff, coupling_mat, res, ~, ~, energy,Cor] = calc_coefficients_nonlinear(cfg);
```
These demo scripts are intended to:

- check that your environment is correctly configured,
- illustrate the expected input data format,
- produce a set of coupling estimates and figures.

## License 
The software in this repository is released under the **Apache License 2.0**, an Open Source Initiative–approved license that permits commercial and non‑commercial use, modification, and redistribution, subject to the terms of the license.

## Citation
If you use this code in your research, please cite the associated manuscript (see the main paper for full citation details) and consider including a link to this repository:

> Schrödinger-Brain: analysis code for complex brain network dynamics.  
> GitHub: https://github.com/ShirleyZhang111/Schrodinger-Brain

## Contact

For questions about the code, bug reports, or suggestions:
- Please open an issue on the GitHub issue tracker of this repository, or  
- Contact the corresponding author as listed in the manuscript.
