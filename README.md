
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bpnreg

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bpnreg)](https://cran.r-project.org/package=bpnreg)
[![Build
Status](https://travis-ci.org/joliencremers/bpnreg.svg?branch=master)](https://travis-ci.org/joliencremers/bpnreg)

The goal of bpnreg is to fit Bayesian projected normal regression models
for circular data.

## Installation

The R-package bpnreg can be installed from CRAN as follows:

``` r
install.packages("bpnreg")
```

You can install a beta-version of bpnreg from github with:

``` r
# install.packages("devtools")
devtools::install_github("joliencremers/bpnreg")
```

## Citation

To cite the package ‘bpnreg’ in publications use:

Jolien Cremers (2020). bpnreg: Bayesian Projected Normal Regression
Models for Circular Data. R package version 1.0.3.
<https://CRAN.R-project.org/package=bpnreg>

## Example

This is a basic example which shows you how to run a Bayesian projected
normal regression model:

``` r
library(bpnreg)
#> Warning: package 'bpnreg' was built under R version 4.0.2
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
#> lppd    -57.02276   8.000000
#> DIC     129.91767   7.933199
#> DIC.alt 129.16896   7.558843
#> WAIC    129.92344   7.938965
#> WAIC2   131.73043   8.842460
#> 
#> 
#> Linear Coefficients 
#> 
#> Component I: 
#>                     mean        mode         sd      LB HPD     UB HPD
#> (Intercept)   1.35903284  1.31054078 0.45916211  0.51176151 2.26408763
#> Condsemi.imp -0.51431167 -0.38356351 0.65112849 -1.76181870 0.77224926
#> Condimp      -0.63880458 -0.74159122 0.66793394 -1.86754185 0.71109682
#> AvAmp        -0.01055016 -0.01139835 0.01218486 -0.03623167 0.01170322
#> 
#> Component II: 
#>                     mean        mode         sd      LB HPD      UB HPD
#> (Intercept)   1.42272991  1.27049889 0.42518913  0.59913060  2.23085396
#> Condsemi.imp -1.17555420 -1.07575082 0.58198521 -2.31181884 -0.04772718
#> Condimp      -0.97477439 -1.16513093 0.61236345 -2.16960668  0.15627423
#> AvAmp        -0.01120924 -0.01173855 0.01088949 -0.03060563  0.01049163
#> 
#> 
#> Circular Coefficients 
#> 
#> Continuous variables: 
#>    mean ax    mode ax      sd ax      LB ax      UB ax 
#>   81.92119   66.77753  106.93783 -126.42475  287.94441 
#> 
#>    mean ac    mode ac      sd ac      LB ac      UB ac 
#>  0.9607303  2.1672257  1.2672149 -0.8206028  2.4573048 
#> 
#>      mean bc      mode bc        sd bc        LB bc        UB bc 
#> -0.001954087  0.011505774  0.030594802 -0.038497492  0.025762454 
#> 
#>      mean AS      mode AS        sd AS        LB AS        UB AS 
#> -0.017532826 -0.007563145  0.309320733 -0.124976392  0.130372158 
#> 
#>    mean SAM    mode SAM      sd SAM      LB SAM      UB SAM 
#> -0.03039514 -0.00759563  0.38996867 -0.25949211  0.22435792 
#> 
#>   mean SSDO   mode SSDO     sd SSDO     LB SSSO     UB SSDO 
#> -0.07185253 -2.05247080  2.04637805 -2.81847543  2.77529945 
#> 
#> Categorical variables: 
#> 
#> Means: 
#>                           mean       mode        sd         LB       UB
#> (Intercept)          0.8133872  0.7973453 0.1984700  0.4171888 1.188380
#> Condsemi.imp         0.2843519  0.2298160 0.4338476 -0.5892295 1.117094
#> Condimp              0.5580093  0.5334042 0.4576417 -0.4067679 1.387852
#> Condsemi.impCondimp -1.2346729 -1.0865211 1.0882791  3.0123533 1.194460
#> 
#> Differences: 
#>                          mean      mode        sd         LB       UB
#> Condsemi.imp        0.5327909 0.5016397 0.5005797 -0.5136915 1.425091
#> Condimp             0.2591943 0.0723281 0.5398564 -0.7832167 1.294051
#> Condsemi.impCondimp 2.1411008 2.4533343 1.0282856 -0.4460970 3.959229
```
