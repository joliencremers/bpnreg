
<!-- README.md is generated from README.Rmd. Please edit that file -->
bpnreg
======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bpnreg)](https://cran.r-project.org/package=bpnreg)
[![Build
Status](https://travis-ci.org/joliencremers/bpnreg.svg?branch=master)](https://travis-ci.org/joliencremers/bpnreg)

The goal of bpnreg is to fit Bayesian projected normal regression models
for circular data.

Installation
------------

You can install bpnreg from github with:

``` r
# install.packages("devtools")
devtools::install_github("joliencremers/bpnreg")
```

Example
-------

This is a basic example which shows you how to run a Bayesian projected
normal regression model:

``` r
library(bpnreg)
bpnr(Phaserad ~ Cond + AvAmp, Motor)
#> Projected Normal Regression 
#> 
#> Model 
#> 
#> Call: 
#> bpnr(pred.I = Phaserad ~ Cond + AvAmp, data = Motor)
#> 
#> MCMC: 
#> iterations = 1000
#> burn-in = 1
#> lag = 1
#> 
#> Model Fit: 
#>         Statistic Parameters
#> lppd    -56.98665   8.000000
#> DIC     130.03264   7.978867
#> DIC.alt 129.00526   7.465178
#> WAIC    130.13423   8.080461
#> WAIC2   131.93160   8.979148
#> 
#> 
#> Linear Coefficients 
#> 
#> Component I: 
#>                     mean         mode         sd      LB HPD     UB HPD
#> (Intercept)   1.38838877  1.430450408 0.44851825  0.54791575 2.21283384
#> Condsemi.imp -0.55387711 -0.586686873 0.61704234 -1.66603600 0.62531778
#> Condimp      -0.64634612 -0.671047696 0.67977534 -1.86378099 0.74237047
#> AvAmp        -0.01081638 -0.007612952 0.01192791 -0.03374693 0.01254156
#> 
#> Component II: 
#>                     mean        mode         sd      LB HPD      UB HPD
#> (Intercept)   1.43186794  1.34887463 0.42859821  0.61193193  2.26798239
#> Condsemi.imp -1.21413507 -1.31468438 0.58965151 -2.29088651 -0.01454310
#> Condimp      -0.97439306 -1.21569705 0.63152408 -2.20567414  0.22262240
#> AvAmp        -0.01174821 -0.01165855 0.01121201 -0.03183777  0.01192664
#> 
#> 
#> Circular Coefficients 
#> 
#> Continuous variables: 
#>    mean ax    mode ax      sd ax      LB ax      UB ax 
#>   91.20303   70.15309  131.17655 -104.32763  313.05562 
#> 
#>    mean ac    mode ac      sd ac      LB ac      UB ac 
#>  0.8709000  2.0828734  1.3282028 -0.8151373  2.5196304 
#> 
#>      mean bc      mode bc        sd bc        LB bc        UB bc 
#> -0.004133268  0.009906482  0.037287197 -0.035000622  0.025243901 
#> 
#>      mean AS      mode AS        sd AS        LB AS        UB AS 
#> -0.015613031 -0.006036166  0.323371494 -0.205164040  0.106346969 
#> 
#>     mean SAM     mode SAM       sd SAM       LB SAM       UB SAM 
#>  0.120674972 -0.008383957  3.571831737 -0.192199588  0.240442932 
#> 
#>   mean SSDO   mode SSDO     sd SSDO     LB SSSO     UB SSDO 
#>  0.02983466 -2.03153425  2.07108578 -2.76494006  2.78617593 
#> 
#> Categorical variables: 
#> 
#> Means: 
#>                           mean       mode        sd         LB       UB
#> (Intercept)          0.8036304  0.7090624 0.1980157  0.4295190 1.180707
#> Condsemi.imp         0.2449676  0.3431420 0.4128835 -0.5628898 1.077003
#> Condimp              0.5505015  0.6278918 0.4600773 -0.4725950 1.344139
#> Condsemi.impCondimp -1.2906544 -1.6598716 1.0721323  3.0402344 1.082592
#> 
#> Differences: 
#>                          mean       mode        sd         LB       UB
#> Condsemi.imp        0.5601076 0.50581982 0.4810523 -0.3500191 1.504551
#> Condimp             0.2564314 0.09884444 0.5388666 -0.8239476 1.297076
#> Condsemi.impCondimp 2.1797273 2.58905512 1.0137416 -0.4643209 3.889132
```
