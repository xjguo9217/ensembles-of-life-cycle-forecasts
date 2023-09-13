## Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts

This repository provides data and R scripts used in our paper "Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts". <br />
<br />
### Data:<br />
Google search interests of social networking websites <br />
Dell computer sales (Acimovic et al. 2018)<br />
<br />
### Forecasting methods:<br />
Bayesian ensembles of single model type: <br />
Bass<br />
Gamma/shifted Gompertz (GSG)<br />
Trapezoid (Trap)<br />
Tilted-Gompertz (TiGo)<br />
Tilted-Gomertz updated using exponential smoothing (TiGo-ETS)<br />
Exponential smoothing state space model with Box-Cox transformation, ARMA errors, Trend and Seasonal components (TBATS)<br />
<br />
Bayesian ensembles of multiple model types: <br />
Bass, GSG, Trap and TiGo<br />
Bass, GSG, Trap, TiGo and TiGo-ETS<br />
<br />
Bayesian updating approaches:<br />
Bass model updated using Extended Kalman Filter with Continuous State and Discrete Observations<br />
Tilted-Gompertz model updated using Kalman Filter<br />
Bayesian Functional Regression (This benchmark model used in the paper is not uploaded here. The code file has been shared privately by the authors of Lei et al. (2023).)<br />
<br />
Other approaches:<br />
Quantile Regression Forest<br />
LightGBM<br />
<br />
### Reference:
Acimovic J, Erize F, Hu K, Thomas DJ, Van Mieghem JA. 2019. Product life cycle data set: Raw and cleaned data of weekly orders for personal computers. Manufacturing & Service Operations Management 21(1) 171–176.

<br />
Lei D, Hu H, Geng D, Zhang J, Qi Y, Liu S, Shen ZM. 2023. New product life cycle curve modeling and forecasting with product attributes and promotion: A bayesian functional approach. Production and Operations Management 32(2) 655–673.


