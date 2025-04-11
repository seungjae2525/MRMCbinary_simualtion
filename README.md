# Simulation Code and Results

This repository contains the simulation code and results for the following paper:
**"Evaluating diagnostic accuracy of binary medical tests in multi-reader multi-case study."**

The materials are organized to reproduce the simulation results in Sections 2.3 and 5 of the manuscript.

There are two main folders, as follows:
- `Simulation in Sec 2.3`: Contains simulation code and results related to Section 2.3 of the manuscript.
- `Simulation in Sec 5`: Contains simulation code and results related to Section 5 of the manuscript.

For details on each folder and the files within it, please refer to the descriptions below.

The simulation study was conducted on a Windows 11 x64 (build 26100) system, utilizing 100 CPU cores through parallel processing, 
and took a total of xx hours and yy minutes to complete.



## Simulation in Sec 2.3

This folder contains R scripts and result files used to generate simulation results for Section 2.3 of the manuscript. 
Below is a brief description of each subfolder and file:


### 📁 Folder Contents

| File                                           | Description                                                             |
|------------------------------------------------|-------------------------------------------------------------------------|
| `Simulation in Sec 2.3.R`                      | R code to perform simulation                                            |
| `Simulation in Sec 2.3.RData`                  | Output files (`.RData`) from `Simulation in Sec 2.3.R`                  |
| `Simulation in Sec 2.3.xlsx`                   | Simulation results (for Tables S4 and S5)                               |
| `Simulation in Sec 2.3_n100_b2.pdf`            | Simulation results (for Figure S1)                                      |
| `Simulation in Sec 2.3_n100_b3.pdf`            | Simulation results (for Figure S2)                                      |
| `Simulation in Sec 2.3_n500_b2.pdf`            | Simulation results (for Figure S3)                                      |
| `Simulation in Sec 2.3_n500_b3.pdf`            | Simulation results (for Figure S4)                                      |
| `Summary and plot for simulation in Sec 2.3.R` | R code to output all simulation results                                 |


### 🔍 Additional Notes
- For detailed simulation settings, please refer to Section 2.3 of the manuscript.
- Make sure to specify an appropriate working directory before running the code.





## Simulation in Sec 5

This folder contains R scripts and result files used to generate simulation results for Section 5 of the manuscript. 
Below is a brief description of each subfolder and file:


### 📁 Folder Contents

| File                                                  | Description                                          |
|-------------------------------------------------------|------------------------------------------------------|
| `Simulation in Sec 5.R`                               | R code to perform simulation                         |
| `Simulation in Sec 5.RData`                           | Output files (`.RData`) from `Simulation in Sec 5.R` |
| `Simulation in Sec 5.xlsx`                            | Simulation results (for Tables S6 and S7)            |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n100_b2.pdf` | Simulation results (for Figure 1)                    |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n100_b3.pdf` | Simulation results (for Figure S5)                   |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n500_b2.pdf` | Simulation results (for Figure S6)                   |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n500_b3.pdf` | Simulation results (for Figure S7)                   |
| `Summary and plot for simulation in Sec 5.R`          | R code to output all simulation results              |


### 🔍 Additional Notes
- For detailed simulation settings, please refer to Section 5 of the manuscript.
- Make sure to specify an appropriate working directory before running the code.
- The currently uploaded results are for the only combination reported in the paper.
- To see results for other combinations, see the results after modifying "beta_0_v", "beta_2k_v", and "beta_3k_v" values ​​in lines 42–44 of `Summary and plot for simulation in Sec 5.R`.


