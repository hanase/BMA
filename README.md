# BMA

[![R build status](https://github.com/hanase/BMA/workflows/R-CMD-check/badge.svg)](https://github.com/hanase/BMA/actions?workflow=R-CMD-check)

R package for Bayesian model averaging and variable selection for linear models,
        generalized linear models and survival models (cox
        regression). 
        
# Bayesian Model Averaging

*Written by Chris Volinsky*
 
Bayesian Model Averaging is a technique designed to help account for the uncertainty inherent in the model selection process, something which traditional statistical analysis often neglects. By averaging over many different competing models, BMA incorporates model uncertainty into conclusions about parameters and prediction. BMA has been applied successfully to many statistical model classes including linear regression, generalized linear models, Cox regression models, and discrete graphical models, in all cases improving predictive performance. Details on these applications can be found in the papers below. 

## Resources

### Statistical Literature

* Most publications related to Bayesian Model Averaging can be found on [Adrian Raftery's research website](https://sites.stat.washington.edu/raftery/Research/bma.html).

* [BMA resources by Mark F.J. Steel](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/steel/steel_homepage/bma)

###  Econometrics Literature
Model combination has been discussed extensively in the econometric literature, usually in the context of combining several experts' forecasts. Bates and Granger (1969) is the forerunner, inspiring a flurry of activity in the field in the early 1970's.

* "The combination of forecasts". J.M. Bates and C.W.J. Granger (1969).
    Operations Research Aquarterly, 20, 451-468.
    [A seminal paper on combining forecasts, inspiring a flurry of activity in the field]

* "Bayesian and non-Bayesian Methods for Combining Models and Forecasts with Applications to Forecasting International Growth Rates". C. Min and A. Zellner (1993). J. of Econometrics, 56, 89-118.
  
* "To Combine or not to Combine? Issues of Combining Forecasts". F.C. Palm and A. Zellner (1992) J. of Forecasting, 11, 687-701.

* "Experience with forecasting univariate time series and the combination of forecasts (with discussion)". P. Newbold and C.W.J. Granger (1974). Journal of the Royal Statistical Society A, 137, 131-149. [Reading in front of RSS. Comments from statisticians are interesting, and mostly negative toward the idea of combining models.]