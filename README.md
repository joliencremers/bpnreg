
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
Models for Circular Data. R package version 2.0.1.
<https://CRAN.R-project.org/package=bpnreg>

## Example

This is a basic example which shows you how to run a Bayesian projected
normal regression model:

``` r
library(bpnreg)
bpnr(Phaserad ~ Cond + AvAmp, Motor, its = 100)
#> Iteration:1
#> Iteration:2
#> Iteration:3
#> Iteration:4
#> Iteration:5
#> Iteration:6
#> Iteration:7
#> Iteration:8
#> Iteration:9
#> Iteration:10
#> Iteration:11
#> Iteration:12
#> Iteration:13
#> Iteration:14
#> Iteration:15
#> Iteration:16
#> Iteration:17
#> Iteration:18
#> Iteration:19
#> Iteration:20
#> Iteration:21
#> Iteration:22
#> Iteration:23
#> Iteration:24
#> Iteration:25
#> Iteration:26
#> Iteration:27
#> Iteration:28
#> Iteration:29
#> Iteration:30
#> Iteration:31
#> Iteration:32
#> Iteration:33
#> Iteration:34
#> Iteration:35
#> Iteration:36
#> Iteration:37
#> Iteration:38
#> Iteration:39
#> Iteration:40
#> Iteration:41
#> Iteration:42
#> Iteration:43
#> Iteration:44
#> Iteration:45
#> Iteration:46
#> Iteration:47
#> Iteration:48
#> Iteration:49
#> Iteration:50
#> Iteration:51
#> Iteration:52
#> Iteration:53
#> Iteration:54
#> Iteration:55
#> Iteration:56
#> Iteration:57
#> Iteration:58
#> Iteration:59
#> Iteration:60
#> Iteration:61
#> Iteration:62
#> Iteration:63
#> Iteration:64
#> Iteration:65
#> Iteration:66
#> Iteration:67
#> Iteration:68
#> Iteration:69
#> Iteration:70
#> Iteration:71
#> Iteration:72
#> Iteration:73
#> Iteration:74
#> Iteration:75
#> Iteration:76
#> Iteration:77
#> Iteration:78
#> Iteration:79
#> Iteration:80
#> Iteration:81
#> Iteration:82
#> Iteration:83
#> Iteration:84
#> Iteration:85
#> Iteration:86
#> Iteration:87
#> Iteration:88
#> Iteration:89
#> Iteration:90
#> Iteration:91
#> Iteration:92
#> Iteration:93
#> Iteration:94
#> Iteration:95
#> Iteration:96
#> Iteration:97
#> Iteration:98
#> Iteration:99
#> Iteration:100
#> Projected Normal Regression 
#> 
#> Model 
#> 
#> Call: 
#> bpnr(pred.I = Phaserad ~ Cond + AvAmp, data = Motor, its = 100)
#> 
#> MCMC: 
#> iterations = 100
#> burn-in = 1
#> lag = 
#> 
#> Model Fit: 
#>         Statistic Parameters
#> lppd     -57.1688   8.000000
#> DIC      127.9570   6.915886
#> DIC.alt  124.5182   5.196498
#> WAIC1    127.7447   6.703544
#> WAIC2    129.1263   7.394339
#> 
#> 
#> Linear Coefficients 
#> 
#> Component I: 
#>                     mean        mode          sd      LB HPD      UB HPD
#> (Intercept)   1.35790309  1.53919307 0.391924091  0.65691407 2.057654675
#> Condsemi.imp -0.52983534 -0.41612729 0.530374398 -1.50572773 0.426828296
#> Condimp      -0.68404666 -0.76754183 0.580782922 -1.65486565 0.289774837
#> AvAmp        -0.01179946 -0.01223479 0.009548015 -0.03090843 0.005276706
#> 
#> Component II: 
#>                     mean         mode          sd     LB HPD     UB HPD
#> (Intercept)   1.42614025  1.079492806 0.416421481  0.6984332  2.2183433
#> Condsemi.imp -1.15627523 -1.063931210 0.538037522 -2.2837229 -0.2885647
#> Condimp      -1.01689511 -1.125072141 0.586648246 -1.9668072  0.1881823
#> AvAmp        -0.01046688 -0.009172757 0.009881872 -0.0306683  0.0055209
#> 
#> 
#> Circular Coefficients 
#> 
#> Continuous variables: 
#>   mean ax   mode ax     sd ax     LB ax     UB ax 
#> 102.35258  73.34450  86.63490  24.19556 367.47488 
#> 
#>    mean ac    mode ac      sd ac      LB ac      UB ac 
#>  0.9268703  1.8524139  1.3298789 -0.7441615  2.4409921 
#> 
#>     mean bc     mode bc       sd bc       LB bc       UB bc 
#> -0.16793096  0.02375924  1.29982126 -0.28692522  0.45828966 
#> 
#>       mean AS       mode AS         sd AS         LB AS         UB AS 
#>  4.380087e-04  3.366778e-05  1.555164e-03 -9.855660e-04  5.396278e-03 
#> 
#>     mean SAM     mode SAM       sd SAM       LB SAM       UB SAM 
#> 2.009564e-04 3.131051e-05 3.626970e-04 7.397841e-06 6.529131e-04 
#> 
#>  mean SSDO  mode SSDO    sd SSDO    LB SSSO    UB SSDO 
#> -0.1083323  1.7910062  2.0399111 -2.8212582  2.5798523 
#> 
#> Categorical variables: 
#> 
#> Means: 
#>                           mean       mode        sd         LB        UB
#> (Intercept)          0.8067426  0.8972646 0.1975172  0.4065758 1.1637551
#> Condsemi.imp         0.2985994  0.1569926 0.3678727 -0.4165081 0.9970036
#> Condimp              0.5623415  0.7778834 0.4861090 -0.4705304 1.3894279
#> Condsemi.impCondimp -1.4038001 -0.9012296 1.1367688  2.5048970 0.8284608
#> 
#> Differences: 
#>                          mean       mode        sd         LB       UB
#> Condsemi.imp        0.5095912  0.3943821 0.4515864 -0.3455296 1.390026
#> Condimp             0.2472478 -0.1522208 0.5688090 -0.9860141 1.138581
#> Condsemi.impCondimp 2.3183579  2.0576422 1.0578694 -0.1311784 4.307274
```
