
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bpnreg

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bpnreg)](https://cran.r-project.org/package=bpnreg)

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

Jolien Cremers (2021). bpnreg: Bayesian Projected Normal Regression
Models for Circular Data. R package version 2.0.2.
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
#> lag = 1
#> 
#> Model Fit: 
#>         Statistic Parameters
#> lppd    -57.22945   8.000000
#> DIC     127.66465   6.768024
#> DIC.alt 124.17298   5.022188
#> WAIC1   127.33436   6.437733
#> WAIC2   128.65389   7.097498
#> 
#> 
#> Linear Coefficients 
#> 
#> Component I: 
#>                      mean        mode         sd      LB HPD     UB HPD
#> (Intercept)   1.319611894  1.39128370 0.45635201  0.33485506 2.03794238
#> Condsemi.imp -0.522451171 -0.47667290 0.57057933 -1.55833243 0.50939839
#> Condimp      -0.650053029 -0.99688228 0.64741848 -2.00362696 0.53197461
#> AvAmp        -0.009320081 -0.01808984 0.01296947 -0.03096035 0.01524266
#> 
#> Component II: 
#>                     mean         mode         sd     LB HPD       UB HPD
#> (Intercept)   1.37081341  1.057909990 0.43448499  0.5256653  2.265534446
#> Condsemi.imp -1.13529041 -1.508829276 0.60583443 -2.2586284  0.029840305
#> Condimp      -0.93550260 -1.263941265 0.62075876 -2.3158274 -0.009041090
#> AvAmp        -0.01016616 -0.003931414 0.01062028 -0.0285245  0.008526117
#> 
#> 
#> Circular Coefficients 
#> 
#> Continuous variables: 
#>    mean ax    mode ax      sd ax      LB ax      UB ax 
#>  116.31973   76.25854  562.60196 -154.19115  219.74298 
#> 
#>    mean ac    mode ac      sd ac      LB ac      UB ac 
#>  1.0746179  2.2543777  1.1994513 -0.8224601  2.4169745 
#> 
#>      mean bc      mode bc        sd bc        LB bc        UB bc 
#> -0.034814814 -0.006854753  0.499046459 -0.767238134  0.666230333 
#> 
#>       mean AS       mode AS         sd AS         LB AS         UB AS 
#>  4.875002e-04  6.466495e-05  5.442953e-03 -1.160784e-02  2.842468e-03 
#> 
#>     mean SAM     mode SAM       sd SAM       LB SAM       UB SAM 
#> 1.437848e-03 1.305745e-04 1.940441e-02 3.180594e-08 3.466995e-03 
#> 
#>   mean SSDO   mode SSDO     sd SSDO     LB SSSO     UB SSDO 
#> -0.05101017  1.88339563  1.99577431 -2.77725635  2.64369230 
#> 
#> Categorical variables: 
#> 
#> Means: 
#>                           mean       mode        sd         LB        UB
#> (Intercept)          0.8119255  0.8675846 0.1957991  0.4326112 1.2082844
#> Condsemi.imp         0.2962062  0.3373583 0.3399843 -0.4996824 0.8360214
#> Condimp              0.5851581  0.4454521 0.4819606 -0.4032866 1.4047517
#> Condsemi.impCondimp -1.3273542 -2.0443304 1.1135480 -2.8870086 1.4407720
#> 
#> Differences: 
#>                          mean      mode        sd         LB       UB
#> Condsemi.imp        0.5152442 0.4826193 0.4033441 -0.2197928 1.286760
#> Condimp             0.2261741 0.3480214 0.5484078 -0.8033373 1.395936
#> Condsemi.impCondimp 2.2043432 2.8593855 1.0362019 -0.4035095 3.855837
```
