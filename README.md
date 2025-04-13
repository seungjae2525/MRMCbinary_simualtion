# Simulation Code and Results

This repository contains the simulation code and results for the following paper:
**"Evaluating diagnostic accuracy of binary medical tests in multi-reader multi-case study."**

The materials are organized to reproduce the simulation results in Sections 2.3 and 5 of the manuscript.

There are two main folders in this repository, as follows:
- `Simulation in Sec 2.3`: Contains simulation code and results related to Section 2.3 of the manuscript.
- `Simulation in Sec 5`: Contains simulation code and results related to Section 5 of the manuscript.

For details on each folder and the files within it, please refer to the descriptions below.



## Simulation in Sec 2.3

This folder contains R scripts, RData, and result (Table and Figure) files used to reproduce simulation results for Section 2.3 in the main manuscript. 

Below is a brief description of each file in this folder:


### 📁 Folder Contents

| File                                           | Description                                            |
|------------------------------------------------|--------------------------------------------------------|
| `Simulation in Sec 2.3.R`                      | R code to perform simulation                           |
| `Simulation in Sec 2.3.RData`                  | Output files (`.RData`) from `Simulation in Sec 2.3.R` |
| `Simulation in Sec 2.3.xlsx`                   | Tables S4 and S5 in the Supplementary Material         |
| `Simulation in Sec 2.3_n100_b2.pdf`            | Figure S1 in the Supplementary Material                |
| `Simulation in Sec 2.3_n100_b3.pdf`            | Figure S2 in the Supplementary Material                |
| `Simulation in Sec 2.3_n500_b2.pdf`            | Figure S3 in the Supplementary Material                |
| `Simulation in Sec 2.3_n500_b3.pdf`            | Figure S4 in the Supplementary Material                |
| `Summary and plot for simulation in Sec 2.3.R` | R code to output all simulation results                |


### 🔍 Additional Notes
- For detailed simulation settings, please refer to Section S2.2 of the Supplementary Material.
- Make sure to specify an appropriate working directory before running the code.
- All simulation studies were conducted on the x86_64-pc-linux-gnu (64-bit) platform, running under Ubuntu 22.04.3 LTS, utilizing 100 CPU cores through parallel processing, and took approximately 13 minutes to complete (for a total of 6 combinations).



## Simulation in Sec 5

This folder contains R scripts, RData, and result (Table and Figure) files used to reproduce simulation results for Section 5 of the manuscript. 

Below is a brief description of each file in this folder:


### 📁 Folder Contents

| File                                                  | Description                                          |
|-------------------------------------------------------|------------------------------------------------------|
| `Simulation in Sec 5.R`                               | R code to perform simulation                         |
| `Simulation in Sec 5.xlsx`                            | Tables S6 and S7 in the Supplementary Material       |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n100_b2.pdf` | Figure 1 in the main manuscript                      |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n100_b3.pdf` | Figure S5 in the Supplementary Material              |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n500_b2.pdf` | Figure S6 in the Supplementary Material              |
| `Simulation in Sec 5_b0(2.5)_b2(1)_b3(1)_n500_b3.pdf` | Figure S7 in the Supplementary Material              |
| `Summary and plot for simulation in Sec 5.R`          | R code to output all simulation results              |

### 📥 Download RData file

- The `Simulation in Sec 5.RData` file size is too large to be uploaded directly to this repository (about 895MB).
- Therefore, you can download the `Simulation in Sec 5.RData` file (output files (`.RData`) from Simulation in Sec 5.R) from the [latest release](https://github.com/seungjae2525/MRMCbinary_simualtion/releases/latest).
- Before loading the downloaded file in R, change the file name from "Simulation.in.Sec.5.RData" to "Simulation in Sec 5.RData".



### 🔍 Additional Notes
- For detailed simulation settings, please refer to Section 5 in the main manuscript.
- Make sure to specify an appropriate working directory before running the code.
- The currently uploaded results are for the only combination reported in the paper.
- To see results for other combinations, see the results after modifying "beta_0_v", "beta_2k_v", and "beta_3k_v" values ​​in lines 42–44 of `Summary and plot for simulation in Sec 5.R`.
- All simulation studies were conducted on the x86_64-pc-linux-gnu (64-bit) platform, running under Ubuntu 22.04.3 LTS, utilizing 100 CPU cores through parallel processing, and took approximately 12 hours to complete (for a total of 288 combinations).


