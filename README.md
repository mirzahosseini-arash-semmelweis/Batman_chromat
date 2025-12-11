# Batman Peaks Analysis with R

This repository contains R scripts and functions for analyzing *Batman peaks* in chromatographic data and determining enantiomerization parameters.  
The workflow implements both the **unified equation** and a **stochastic model**, providing a flexible and automated approach for extracting kinetic and thermodynamic constants from chromatograms.

---

## Features

- Automated **peak picking** using the [`pracma`](https://cran.r-project.org/package=pracma) pattern recognition algorithm.
- Support for both **unified equation** and **stochastic modeling** (using the [`DEoptim`](https://cran.r-project.org/web/packages/DEoptim/index.html) and [`minpack.lm`](https://cran.r-project.org/web/packages/minpack.lm/index.html) packages) of Batman peaks.
- **Baseline correction** with polynomial fitting.
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

---

## Tips on Use

* Place all chromatogram files (`.CSV`, containing **two columns: time, intensity**) into a single folder.
* Name the .csv files with the following convention for the algorithm to fetch correct flow rate and temperature parameters: **`[compound]_[column]_[flow]_[temperature].csv`**
* Run **`plot.runs()`** to visually inspect the chromatograms and verify signal quality. Note that the algorithm assumes that if void time peak is present, it is discernible (free from artefacts; required for deconvolution from chromatogram) and it has lower intensity than the Batman peak (required for automatic peak detection).
* Execute **`batch.eval.kue()`** with default parameters to perform automated unified equation fitting, and inspect the generated plots for consistency.
* If some peaks are not detected or misidentified:
  * Adjust **`threshold`** (broader peaks often require a lower threshold, maybe even **`0`**).
  * Modify **`minSNR`** so that even the lowest-SNR peak remains identifiable.
  * In rare cases where peaks are very close, set **`minpeakdist = 0`**.
* If peak identification still needs refinement, create a CSV file named **`peak_table.csv`** with the following columns:
  ```
  file_name, t_M, t_A, t_B
  ```
  (Use file names exactly as listed in the output `summary_data.csv`.)
  Provided *t* values should correspond to the expected peak positions.
* Re-run **`batch.eval.kue()`** with the argument **`peak_table = TRUE`** to use your manually verified peak positions. Adjust values in `peak_table.csv` if necessary until all runs are correctly identified.
* For **coalesced (merged) peaks** that cannot be evaluated using the unified equation:
  * Leave the columns **`t_A`** and **`t_B`** empty in `peak_table.csv`.
  * This signals the algorithm to apply only the **stochastic model** to those runs.
  * However, providing **`t_M`** (the void time) is still recommended for accurate stochastic evaluation.
* (Optional) If the **non-retained analyte** (used for determining the void time) is **not present in every run**, but measured separately:
  * Provide these values in **`peak_table.csv`** as **`t_M`**.
  * Then, run **`batch.eval.kue()`** with **`t0_given = TRUE`**.
* Finally, perform the **stochastic model fitting** by running:

  ```r
  batch.eval.stoch()
  ```

## Script Versions

* `Batman_script.R` contains code published in [Journal of Chromatography A](https://doi.org/10.1016/j.chroma.2025.466558), originally implemented unified equation and stochastic model fit with one-site model of chiral stationary phase
* `Batman_v2_script.R` contains code with implementation of two-site model of chiral stationary phase in the stochastic model as well. By setting **`two_site = TRUE`** the algorithm will fit the extended model (otherwise defaults to one-site model), if nonselective/enantioselective site ratio in log scale is provided as dataframe for all *temperature × stationary phase* combinations present in the experiment (**`log_ratio`**).
