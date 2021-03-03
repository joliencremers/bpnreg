# bpnreg 2.0.0.9000

* This is a development version

# bpnreg 2.0.0
## Major changes
* Changed the names of the mcmc output objects. `B1`, `B2`, `Beta.I`, `Beta.II`, `B.I`, `B.II`, `VCov.I` and `VCov.II` are now named `beta1`, `beta2`, `beta1`, `beta2`, `b1`, `b2`, `omega1` and `omega2` respectively.

* Functions `predict()` and `residual()` were removed.

* MCMC sampler underlying the `bpnme` function is rewritten in `c++`. This should significantly increase computational speed.

* Increased memory efficiency by not storing intermediate objects in memory.

* A vignette containing answers to frequently asked questions about usage of the package is included.

## Bug fixes
* Adapted documentation `bpnme` function such that it clearly states only one grouping factor can be used.

* Made sure that the circular random intercept variance (and circular random slope variance of categorical predictors) is computed as 1 - mean resultant length. Also updated documentation accordingly.

* Results for `bc`, `AS` and `SAM`  are now also included in degrees if `units = "degrees"` in `circ_coef` function

* An error warning in case of any missings in the data was included.

* The `WAIC` criterion in the model fit estimates obtained from the `fit()` function is more correctly referred to as `WAIC1`.

* Format of `cRS` object stored in a `bpnme` fit object was changed from `matrix` to `data.frame` to accomodate for data of more than two types, e.g. `character` and `numeric`.

# bpnreg 1.0.3
* Deleted arguments from function documentation that are not in `\usage`.

* Updated references in the packages description.

# bpnreg 1.0.2
* This version contains an update of code containing the class function to make sure the package works with R version 4.0.0

* The check for the range of the circular outcome was corrected.

# bpnreg 1.0.1

* The function bpnme can now also handle the standard R data.frame instead of only Data Frame Tbl (from dplyr, tbl_df(data))

* The function coef.bpnr now correctly returns radians instead of degrees when option units == "radians" is chosen.

* The documentation now states that the circular outcome needs to be measured in radians.

* The references in the documentation have been updated.

* The functions mmr and mmme now give out a warning message when the circular outcome is out of range.

* An error in the assignment of column names in the random effect model matrix (in mmme function) was fixed.

* An error in the check for whether random effects are numeric (in bpnme function) was fixed.

* Several checks and error warnings for wrong specifications of nesting structure were added.

* Internal functions were taken out of the package documentation on CRAN.

* The priors for the (fixed effect) regression coefficients in the regression and mixed-effects model are now both N(0, 10000), We also note this in the documentation of bpnr() and bpnme().

* We note the fact that the mixed-effects model is only developed for models with one nesting variable in the documentation.

# bpnreg 1.0.0
* This is the first version of the package.
