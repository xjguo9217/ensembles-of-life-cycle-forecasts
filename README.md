## Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts

This repository provides data and R scripts used in our paper "Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts". Below is the list and descriptions of the files.

### Study 1: Social Networks 


### Study 2: Computer Sales

| Data set| Descriptions |
| --- | --- | 
| dell_data_raw.csv  | Normalized sales without end-of-life orders truncated (Acimovic et al. 2018). This data is used to simulate net cumulative sales (M). |
| dell_data_truncated.csv| Standardized sales with end-of-life orders truncated. The truncated time series are used in Hu et al.(2018) for forecasting computer sales.|
| weekly_M.csv | A data file created by running "S1 - SampleCumulativeSales - M.R". This file stores the simulated net cumulative sales (M).


#### Code files
| Code file | Descriptions |
| --- | --- | 
| S1 - SampleCumulativeSales - M.R  | Data preparation in Study 2. Sample net cumulative sales (M) based on the method described in Section A.7.|




