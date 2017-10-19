
<!-- README.md is generated from README.Rmd. Please edit that file -->
bpnreg
======

The goal of bpnreg is to fit Bayesian projected normal regression models for circular data.

Installation
------------

You can install bpnreg from github with:

``` r
# install.packages("devtools")
devtools::install_github("joliencremers/bpnreg")
```

Example
-------

This is a basic example which shows you how to run a Bayesian projected normal regression model:

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
#> lppd     -57.1423   8.000000
#> DIC      130.1116   8.024107
#> DIC.alt  132.1696   9.053110
#> WAIC     129.8904   7.802881
#> WAIC2    131.6491   8.682262
#> 
#> 
#> Linear Coefficients 
#> 
#> Component I: 
#>                      mean         mode         sd      LB HPD     UB HPD
#> (Intercept)   1.281271224  1.434526732 0.45252272  0.46317178 2.23113550
#> Condsemi.imp -0.408424353 -0.441216754 0.62526234 -1.56249571 0.79919832
#> Condexp      -0.518877035 -0.524355136 0.64099499 -1.74163992 0.73733390
#> AvAmp        -0.008747752 -0.007506795 0.01135152 -0.03100102 0.01353229
#> 
#> Component II: 
#>                     mean        mode         sd      LB HPD      UB HPD
#> (Intercept)   1.39496052  1.49923392 0.43825607  0.54914605 2.199253307
#> Condsemi.imp -1.13706805 -1.01473377 0.60483330 -2.28899755 0.011343999
#> Condexp      -0.97565953 -0.82912074 0.62703430 -2.24519673 0.218029789
#> AvAmp        -0.01144042 -0.01193592 0.01103678 -0.03361457 0.009573006
#> 
#> 
#> Circular Coefficients 
#> 
#> Continuous variables: 
#>    mean ax    mode ax      sd ax      LB ax      UB ax 
#>   81.60681   74.10766  115.26201 -154.13426  276.10493 
#> 
#>    mean ac    mode ac      sd ac      LB ac      UB ac 
#>  0.6509619 -0.4380238  1.2648523 -0.7926655  2.5384750 
#> 
#>     mean bc     mode bc       sd bc       LB bc       UB bc 
#> -0.01364625  0.01066514  0.19358642 -0.04545137  0.02442152 
#> 
#>      mean AS      mode AS        sd AS        LB AS        UB AS 
#> -0.003992539 -0.006218847  0.458896772 -0.117588533  0.164513332 
#> 
#>     mean SAM     mode SAM       sd SAM       LB SAM       UB SAM 
#> -0.012349538 -0.006010997  0.392240548 -0.218694004  0.145853579 
#> 
#>  mean SSDO  mode SSDO    sd SSDO    LB SSSO    UB SSDO 
#>  0.2523587  1.7567794  1.9639745 -2.6364423  2.8579686 
#> 
#> Categorical variables: 
#>              mean circDiff mode circDiff sd circDiff LB circDiff
#> Condsemi.imp     0.5539928     0.6791023   0.4694104  -0.3482197
#> Condexp          0.3397840     0.3974468   0.5334749  -0.7485248
#>              UB circDiff
#> Condsemi.imp    1.516208
#> Condexp         1.403599
```
