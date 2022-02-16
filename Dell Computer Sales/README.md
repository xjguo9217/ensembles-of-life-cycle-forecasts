## Forecasting Dell Computer Sales

This folder contains the following R scripts and data files:

GoogleTrendsTimeSeriesSocial-May-2021.csv: Search interests of 161 social networking websites. <br />
WikiListSocial-startenddays-May-2021.csv: Start and end dates of the social networking websites. <br />
<br /> 
Source code: This folder contains the source code for fitting and predicting from the Bass, GSG, Trap and TiGo-ETS models.<br />
<br /> 
Social networks - Ensemble of diffusion models: Produce Bayesian ensemble forecasts using one of the following model types: Bass, GSG, Trap or TiGo (time-invariant).<br /> 
Social networks - Ensemble of TiGo-ETS: Produce Bayesian ensemble forecasts using the TiGo-ETS model. This is the time-varying version of the TiGo model. <br /> 
Social networks - Bass-EKF: Produce forecasts using the Bass model updated by Extended Kalman Filter with Continuous State and Discrete Observations. This is the time-varying version of the Bass model. <br /> 
Social networks - TiGo-KF: Produce forecasts using the TiGo model updated by Kalman Filter.<br /> 
Social networks - Machine learning models: Produce forecasts using the Quantile Regression Forest or LightGBM.<br /> 
Social networks - Ensemble of four TI diffusion models: Produce Bayesian ensemble forecasts using multiple model types (Bass, GSG, Trap and TiGo (time-invariant))<br /> 
Social networks - Ensemble of all five models: Produce Bayesian ensemble forecasts using multiple model types (Bass, GSG, Trap, TiGo (time-invariant) and TiGo-ETS)<br /> 



