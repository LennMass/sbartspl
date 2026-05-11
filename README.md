# SBART+SPL for PATE

This repository contains replicaton code for the project:

**"Estimating PATE under positivity violations: SBART+SPL for
high-dimensional covariates"**

CLeaR 2026, forthcoming in the [Proceedings of Machine Learning Research 323](https://openreview.net/pdf?id=v6gn5xdbu3).

---

## Project Structure

```
sbartspl/
|
|-- 00_config.R
|
|-- 00_data/
|   |-- 01_simulations/
|       |-- 01_main/
|       |-- 02_additional/
|           |-- 01_dgp_bcf/
|           |-- 02_varinf_param/
|   |-- 02_applications/
|
|-- 01_code/
|   |-- 00_functions/
|       |-- 01_helper/
|       |-- 02_methods/
|   |-- 01_simulations/
|       |-- 01_main/
|       |-- 02_additional/
|           |-- 01_dgp_bcf/
|           |-- 02_varinf_param/
|   |-- 02_application/
|       |-- natgas_leukemia/
|       |-- natgas_thycancer/
|       |-- visual_counties/
|   |-- 03_motivation/
|
|-- 02_results/
|   |-- 01_simulations/
|       |-- 01_main/
|       |-- 02_additional/
|           |-- 01_dgp_bcf/
|           |-- 02_varinf_param/
|   |-- 02_application/
|       |-- natgas_leukemia/
|       |-- natgas_thycancer/
|       |-- visual_counties/
|
|-- .gitignore
|-- sbartspl.Rproj
```


---

## Folder and File Description

### 0_config.R

- Project-wide configuration file (loads packages and self-written functions)


```r
source("0_config.R")
```


---

### 00_data/

- Raw data for main simulations based on Nethery et al. (2019), the subsequent simulations based on Wang et al. (2024), and illustrations for the variance inflation parameter 

---

### 01_code/

- Helper functions (i.e., Bayesian bootstrap, DGPs, region of overlap definition)
- Functions that implement methods (i.e., SBART+SPL, BART+SPL of Nethery et al. (2019), GR of Gutman and Rubin (2013, 2015))
- Data generation, effect estimation and reporting files for simulations and illustrations
- Files to replicate empirical application


---

### 02_results/


- Results for the simulations and empirical application

---




