## Forecasting Dell Computer Sales

This folder contains the following R scripts and data files:

dell_data_raw.csv: Raw sales data of Dell computers. This data is used for sampling cumulative sales M. <br />
dell_data_truncated.csv: Truncated end-of-life sales of Dell computers. This data is used as the life cycles of the computers.
The above two data sets were published in Acimovic et al. (2018).<br />
weekly_M.csv: Sampled cumulative sales M from running the R file "Computer sales - Sample cumulative sales M.R". <br />
<br /> 
Source code: This folder contains the source code for fitting and predicting from the Bass, GSG, Trap and TiGo-ETS models.<br />
<br /> 
Computer sales - Sample cumulative sales M.R: Data Preparation for the Study of Computer Sales. Sample M to unnormalize the normalized sales data. See Section A.7 in the paper for more detail.<br /> 
Computer sales - Ensemble of diffusion models: Produce Bayesian ensemble forecasts using one of the following model types: Bass, GSG, Trap or TiGo (time-invariant).<br /> 
Computer sales - Ensemble of TiGo-ETS: Produce Bayesian ensemble forecasts using the TiGo-ETS model. This is the time-varying version of the TiGo model. <br /> 
Computer sales - Bass-EKF: Produce forecasts using the Bass model updated by Extended Kalman Filter with Continuous State and Discrete Observations. This is the time-varying version of the Bass model. <br /> 
Computer sales - TiGo-KF: Produce forecasts using the TiGo model updated by Kalman Filter.<br /> 
Computer sales - Machine learning models: Produce forecasts using the Quantile Regression Forest or LightGBM.<br /> 
Computer sales - Ensemble of four TI diffusion models: Produce Bayesian ensemble forecasts using multiple time-invariant model types (Bass, GSG, Trap and TiGo (time-invariant))<br /> 
Computer sales - Ensemble of all five models: Produce Bayesian ensemble forecasts using multiple model types (Bass, GSG, Trap, TiGo (time-invariant) and TiGo-ETS)<br /> 

### Reference
Acimovic J, Erize F, Hu K, Thomas DJ, Van Mieghem JA. 2018. Product life cycle data set: Raw andcleaned data of weekly orders for personal computers.Manufacturing & Service Operations Management
