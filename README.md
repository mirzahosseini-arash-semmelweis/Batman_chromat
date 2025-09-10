# Batman Peaks Analysis with R

This repository contains R scripts and functions for analyzing *Batman peaks* in chromatographic data and determining enantiomerization kinetics.  
The workflow implements both the **unified equation** and a **stochastic model**, providing a flexible and automated approach for extracting kinetic and thermodynamic constants from chromatograms.

---

## Features

- Automated **peak picking** using the [`pracma`](https://cran.r-project.org/package=pracma) pattern recognition algorithm.  
- Support for both **unified equation** and **stochastic modeling** of Batman peaks.  
- **Baseline correction** and spline/polynomial fitting.  
- Calculation of kinetic rate constants (`kue_f`, `kue_r`) and thermodynamic parameters (ΔH‡, ΔS‡, ΔG‡).  
- Visualization functions (chromatograms, Eyring–Polanyi plots).  
- Fully **customizable R functions** for extending to new data sets.

---

## Requirements

- [R](https://cran.r-project.org/) (>= 4.0.0)  
- [RStudio](https://posit.co/) (recommended)  
- Required R packages:  
  ```r
  install.packages(c("pracma", "minpack.lm", "broom", "DEoptim", "parallel", "gamlss.dist"))
