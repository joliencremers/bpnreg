# bpnreg 1.0.2
This version contains an update of code containing the class function to make sure the package works with R version 4.0.0

The check for the range of the circular outcome was corrected.

# bpnreg 1.0.1

The function bpnme can now also handle the standard R data.frame instead of only Data Frame Tbl (from dplyr, tbl_df(data))

The function coef.bpnr now correctly returns radians instead of degrees when option units == "radians" is chosen.

The documentation now states that the circular outcome needs to be measured in radians.

The references in the documentation have been updated.

The functions mmr and mmme now give out a warning message when the circular outcome is out of range.

An error in the assignment of column names in the random effect model matrix (in mmme function) was fixed.

An error in the check for whether random effects are numeric (in bpnme function) was fixed.

Several checks and error warnings for wrong specifications of nesting structure were added.

Internal functions were taken out of the package documentation on CRAN.

The priors for the (fixed effect) regression coefficients in the regression and mixed-effects model are now both N(0, 10000), We also note this in the documentation of bpnr() and bpnme().

We note the fact that the mixed-effects model is only developed for models with one nesting variable in the documentation.

# bpnreg 1.0.0
This is the first version of the package.
