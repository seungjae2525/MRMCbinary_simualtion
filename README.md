# Simulation Code and Results

This repository contains the simulation code and results for the following paper:
**"Evaluating diagnostic accuracy of binary medical tests in multi-reader multi-case study."**

The materials are organized by the sections of the manuscript that they support.

There are two main folders:
- `Simulation in Sec 2.3`: Contains simulation code and results related to Section 2.3 of the manuscript.
- `Simulation in Sec 5`: Contains simulation code and results related to Section 5 of the manuscript.

For details on each folder, please refer to the descriptions below.




## Simulation in Sec 2.3

This folder contains R scripts and result files used to generate simulation results for Section 2.3 of the manuscript. Below is a brief description of each subfolder and file:

### 📁 Folder Contents

| File                                           | Description                                                             |
|------------------------------------------------|-------------------------------------------------------------------------|
| `Simulation in Sec 2.3.R`                      | Main scripts for simulation: data generation, model fitting, summary.   |
| `Simulation in Sec 2.3.RData`                  | Output files (`.rds`) with estimates, convergence info, timing.         |
| `Simulation in Sec 2.3.xlsx`                   | Overview of folder contents.                                            |
| `Simulation in Sec 2.3_n100_b2.pdf`            | Overview of folder contents.                                            |
| `Simulation in Sec 2.3_n100_b3.pdf`            | Overview of folder contents.                                            |
| `Simulation in Sec 2.3_n500_b2.pdf`            | Overview of folder contents.                                            |
| `Simulation in Sec 2.3_n500_b3.pdf`            | Overview of folder contents.                                            |
| `Summary and plot for simulation in Sec 2.3.R` | Creates summary tables and plots for Section 2.3 simulation.            |



### 🔍 Additional Notes

- Make sure to run the scripts in the order described in `Rcodes/main.R` if you wish to reproduce the simulations.
- Simulation settings (e.g., number of replications, parameter values) are defined in `Rcodes/setting_sim.R`.




## Simulation in Sec 5

This folder contains R scripts and result files used to generate simulation results for Section 2.3 of the manuscript. Below is a brief description of each subfolder and file:

### 📁 Folder Contents

| File | Description |
|---------------|-------------|
| `Rcodes`      | Contains all R scripts used for the simulation study. These include the main loop, data generation, model fitting functions, and summary routines. |
| `Results`     | Output files (e.g., `.rds`) generated from the simulation. These include saved estimates, convergence status, and timing information. |
| `README.md`   | This file. Provides an overview of the contents in this folder. |

### 🔍 Additional Notes

- Make sure to run the scripts in the order described in `Rcodes/main.R` if you wish to reproduce the simulations.
- Simulation settings (e.g., number of replications, parameter values) are defined in `Rcodes/setting_sim.R`.
