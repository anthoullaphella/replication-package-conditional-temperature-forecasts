# Taking the Highway or the Green Road?  
**Conditional Temperature Forecasts Under Alternative SSP Scenarios: Replication Package**  

**Authors:** Anthoulla Phella∗†, Vasco J. Gabriel‡, Luis F. Martins§  
**Date:** January 2026  

∗Corresponding author: anthoulla.phella@glasgow.ac.uk  
†University of Glasgow  
‡University of Victoria and NIPE  
§ISCTE - Instituto Universitário de Lisboa and CIMS  

---

## 1. Overview

This replication package contains the data and MATLAB code required to reproduce the empirical results, figures, and forecasting exercises presented in the paper.  

The package allows replication of:

- Conditional and unconditional temperature forecasts under alternative SSP scenarios  
- Counterfactual forecasting exercises  
- Forecast evaluation results  

---

## 2. Data Availability and Description

The package includes the following datasets:

- `Final Dataset Forecast SSP.xlsx`  
  Contains temperature anomalies and climate forcings data; observed data until 2023 and SSP projections up to 2050. Used for real-time forecasting.

- `Final Dataset Forecast SSP Counterfactual.xlsx`  
  Contains temperature anomalies and climate forcings data; observed data until 2016 and SSP projections up to 2050. Used for counterfactual forecasting.

- `Out-of-Sample Forecasts.xlsx`  
  Contains results of the forecasting evaluation exercise, including model comparisons.  

- `Forecasts Counterfactual.mat` (intermediary dataset)  
  Contains conditional and unconditional forecasts used in the pseudo out-of-sample evaluation.

All data required to reproduce the results are included in this package.  

---

## 3. Data Processing

The datasets included in this package are final and ready for use. No additional preprocessing is required. All scripts directly load the provided data files.

---

## 4. Computational Requirements

- **Software:** MATLAB  
- **Version:** MATLAB R2021b  
- **Additional requirements:** Custom functions located in the replication package folder `Functions`

**Expected runtime:**

- `mainForecast.m`: approximately 15 minutes  
- `mainForecastCounterfactual.m`: approximately 15 minutes  
- `ForecastEvaluation.m`: less than a minute  
- `PathwaysPlotting.m`: less than a minute  

> Forecasts from the proposed model are generated using `mainForecastCounterfactual.m`, while alternative model forecasts are obtained from publicly available reports by the World Meteorological Organization (WMO) Lead Centre for Annual-to-Decadal Climate Prediction, hosted by the UK Met Office.

---

## 5. Mapping of Code to Results

| Outputs | Script | Notes / Data Used |
|---------|--------|------------------|
| Figure 1 | PathwaysPlottimg.m | Plots SSP scenario pathways for three conditioning variables. Uses `Final Dataset Forecast SSP.xlsx`. Produces equality and inequality constraints for CO2, CH4, N2O under adverse and optimistic scenarios up to 2050. |
| Figure 2 | mainForecast.m | Real-time forecast of temperature anomalies under adverse scenario. Uses `Final Dataset Forecast SSP.xlsx`. Conditional and unconditional forecasts when CO2, CH4, N2O emissions in 2024–2050 match SSP adverse scenario projections. |
| Figure 3 | mainForecast.m | Difference between conditional and unconditional forecasts for adverse scenario. Uses `Final Dataset Forecast SSP.xlsx`. |
| Figure 4 | mainForecast.m | Conditional and unconditional forecasts under SSP optimistic scenario. Uses `Final Dataset Forecast SSP.xlsx`. |
| Figure 5 | mainForecast.m | Difference between conditional and unconditional forecasts under SSP optimistic scenario. Uses `Final Dataset Forecast SSP.xlsx`. |
| Figure 6 | mainForecastCounterfactual.m | Counterfactual forecast of temperature anomalies under optimistic scenario achieved since 2016. Uses `Final Dataset Forecast SSP Counterfactuals.xlsx`. Conditional and unconditional forecasts for CO2, CH4, N2O emissions in 2016–2050. |
| Figure 7 | mainForecastCounterfactual.m | Difference between conditional and unconditional forecasts for counterfactual scenario. Uses `Final Dataset Forecast SSP Counterfactuals.xlsx`. |
| Figures 8–15 | mainForecast.m | Appendix figures under alternative priors. Uses `Final Dataset Forecast SSP.xlsx`. |
| Table 2 | ForecastEvaluation.m | Runs pseudo out-of-sample forecast evaluation exercise. Uses `Out-of-Sample Forecasts.xlsx`. |

---

## 6. License

The code in this replication package is licensed under the **MIT License**.  

The data are derived from publicly available sources and are subject to their original licenses. Where applicable, processed datasets included in this package are shared under the **Creative Commons Attribution 4.0 (CC BY 4.0)** license.

---

## 7. Notes

- Ensure that all functions in the `/Functions` folder are in your MATLAB path. You can add them using:

```matlab
addpath(genpath('Functions'))
