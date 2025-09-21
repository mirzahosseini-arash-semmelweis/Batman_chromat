# Batman Peaks Analysis with R

This repository contains R scripts and functions for analyzing *Batman peaks* in chromatographic data and determining enantiomerization parameters.  
The workflow implements both the **unified equation** and a **stochastic model**, providing a flexible and automated approach for extracting kinetic and thermodynamic constants from chromatograms.

---

## Features

- Automated **peak picking** using the [`pracma`](https://cran.r-project.org/package=pracma) pattern recognition algorithm.
- Support for both **unified equation** and **stochastic modeling** (using the [`DEoptim`](https://cran.r-project.org/web/packages/DEoptim/index.html) and [`minpack.lm`](https://cran.r-project.org/web/packages/minpack.lm/index.html) packages) of Batman peaks.
- **Baseline correction** and spline/polynomial fitting.
- Calculation of kinetic rate constants (`kue_f`, `kue_r`) and thermodynamic parameters (ΔH‡, ΔS‡, ΔG‡).
- Visualization functions of chromatograms and model fits.
- Fully **customizable R functions** for extending to new data sets.

---

## Requirements

- [R](https://cran.r-project.org/) (>= 4.0.0)  
- [RStudio](https://posit.co/) (recommended)  
- Required R packages:  
  ```r
  install.packages(c("pracma", "minpack.lm", "broom", "DEoptim", "parallel", "gamlss.dist"))

---

## File format

- Chromatograms should be provided as individual .CSV files, each containing two columns: time and intensity, and placed together in a designated folder
- If your chromatograms are in another format, they can be converted to .CSV using the [`chromConverter`](https://cran.r-project.org/web/packages/chromConverter/index.html) package
- For multichannel chromatograms, the user must select the **optimal detector signal** prior to conversion
