## Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts

This repository provides data and R scripts used in our paper "Bayesian Ensembles of Exponentially Smoothed Life-Cycle Forecasts". <br />
<br />
Data:<br />
Google search interests of social networking websites <br />
Dell computer sales (Acimovic et al. 2018)<br />
<br />
Forecasting methods:<br />
Bayesian ensembles of single model type: <br />
1.  Bass<br />
2.  Gamma/shifted Gompertz (GSG)<br />
3.  Trapezoid (Trap)<br />
4.  Tilted-Gompertz (TiGo)<br />
5.  Tilted-Gomertz updated using exponential smoothing (TiGo-ETS)<br />
<br />
Bayesian ensembles of multiple model types: <br />
6.  Bass, GSG, Trap and TiGo<br />
7.  Bass, GSG, Trap, TiGo and TiGo-ETS<br />
<br />
Other approaches:<br />
8.  Bass model updated using Extended Kalman Filter with Continuous State and Discrete Observations<br />
9.  Tilted-Gompertz model updated using Kalman Filter<br />
10. Quantile Regression Forest<br />
11. LightGBM<br />
<br />
### Reference
Acimovic J, Erize F, Hu K, Thomas DJ, Van Mieghem JA. 2018. Product life cycle data set: Raw and cleaned data of weekly orders for personal computers. Manufacturing & Service Operations Management

