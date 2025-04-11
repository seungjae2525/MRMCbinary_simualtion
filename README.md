# MRMCbinary_simualtion

This directory contains the simulation code and results for the paper  
**"Evaluating diagnostic accuracy of binary medical tests in multi-reader multi-case study."**


## Simulation in Section 2.3

This folder contains R scripts and result files used to generate simulation results for Section 2.3 of the manuscript. Below is a brief description of each subfolder and file:

### 📁 Folder Contents

| Simulation in Sec 2.3 / File | Description |
|---------------|-------------|
| `Rcodes`      | Contains all R scripts used for the simulation study. These include the main loop, data generation, model fitting functions, and summary routines. |
| `Results`     | Output files (e.g., `.rds`) generated from the simulation. These include saved estimates, convergence status, and timing information. |
| `README.md`   | This file. Provides an overview of the contents in this folder. |

### 🔍 Additional Notes

- Make sure to run the scripts in the order described in `Rcodes/main.R` if you wish to reproduce the simulations.
- Simulation settings (e.g., number of replications, parameter values) are defined in `Rcodes/setting_sim.R`.




## Simulation in Section 5

This folder contains R scripts and result files used to generate simulation results for Section 2.3 of the manuscript. Below is a brief description of each subfolder and file:

### 📁 Folder Contents

| Simulation in Sec 5 / File | Description |
|---------------|-------------|
| `Rcodes`      | Contains all R scripts used for the simulation study. These include the main loop, data generation, model fitting functions, and summary routines. |
| `Results`     | Output files (e.g., `.rds`) generated from the simulation. These include saved estimates, convergence status, and timing information. |
| `README.md`   | This file. Provides an overview of the contents in this folder. |

### 🔍 Additional Notes

- Make sure to run the scripts in the order described in `Rcodes/main.R` if you wish to reproduce the simulations.
- Simulation settings (e.g., number of replications, parameter values) are defined in `Rcodes/setting_sim.R`.
