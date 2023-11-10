## Forecasting Demand of Intel Microprocessors
Code for the Bayesian Functional Regression used in the paper is not uploaded here. The code file has been shared privately by the authors of Lei et al. (2023).)

This folder contains the following R scripts and data files:

Intel_actual_demand.csv: Raw sales data of Intel microprocessors.  <br />
startend_date.csv: Start and end week of the product. <br />
The above two data sets were published in Manary and Willems (2022).<br />
<br /> 
Source code: This folder contains the source code for fitting and predicting from the Bass, GSG, Trap and TiGo-ETS models.<br />
<br /> 
Intel sales - Ensemble of diffusion models: Produce Bayesian ensemble forecasts using one of the following model types: Bass, GSG, Trap or TiGo (time-invariant).<br /> 
Intel sales - Ensemble of TiGo-ETS: Produce Bayesian ensemble forecasts using the TiGo-ETS model. This is the time-varying version of the TiGo model. <br /> 
Intel sales - Bass-EKF: Produce forecasts using the Bass model updated by Extended Kalman Filter with Continuous State and Discrete Observations. This is the time-varying version of the Bass model. <br /> 
Intel sales - TiGo-KF: Produce forecasts using the TiGo model updated by Kalman Filter.<br /> 
Intel sales - Machine learning models: Produce forecasts using the Quantile Regression Forest or LightGBM.<br /> 
Intel sales - Ensemble of four TI diffusion models: Produce Bayesian ensemble forecasts using multiple time-invariant model types (Bass, GSG, Trap and TiGo (time-invariant))<br /> 
Intel sales - Ensemble of all five models: Produce Bayesian ensemble forecasts using multiple model types (Bass, GSG, Trap, TiGo (time-invariant) and TiGo-ETS)<br /> 
Intel sales - Bayesian nonparametric: Produce forecasts from the Bayesian nonparametric model proposed in Dew et al. (2018)

### Reference
Manary MP, Willems SP. 2022. Data set: 187 weeks of customer forecasts and orders for microprocessors from intel corporation. Manufacturing & Service Operations Management 24(1) 682–689.

Dew R, Ansari A. 2018. Bayesian nonparametric customer base analysis with model-based visualizations. Marketing Science 37(2) 216–235.
