---
affiliations:
- index: 1
  name: |
      Section of Biostatistics, Department of Public Health, University of
      Copenhagen
authors:
- affiliation: 1
  name: Jolien Cremers
  orcid: '0000-0002-0037-1747'
bibliography: 'paper.bib'
date: '24/02/2020'
output:
  md_document:
    pandoc_args: '--atx-headers'
    preserve_yaml: True
    variant: markdown
  pdf_document: default
tags:
- R
- Bayesian statistics
- Circular data
- Regression
- 'Mixed-effects models'
title: |
    bpnreg: A package to analyze Bayesian projected normal circular
    regression models
---

# Summary

# Statement of Need

The purpose of `bpnreg` is to provide methods for fitting circular
regression and mixed-effects models to `R`-users. THe main functions
allow for fitting Bayesian multiple and mixed-effect regression models
for circular data based on the projected normal distribution (see
@Nunez-Antonio2011-fm and @Nunez-Antonio2014-bd for a description of the
models). Both continuous and categorical predictors can be included.
Sampling from the posterior is performed via an MCMC algorithm
implemented in `c++` that allows for fast computation (see
@Cremers2018-ta and @Cremers2021-mm for a description of the MCMC
samplers). Posterior descriptives of all parameters, model fit
statistics and Bayes factors for hypothesis tests for inequality
constrained hypotheses are provided.

Only a couple of R-packages provide methods for circular data analyses
and even less contain functionality for circular regression models. An
overview is given in @pewsey2020recent. The package `circular`
[@Agostinelli2017] is a general purpose package that also contains
functionality for fitting frequentist regression models for circular
outcomes based on the von Mises distribution. The package `circglmbayes`
[@Mulder2017] provides a Bayesian regression model for circular outcomes
based on the von Mises distribution. To date `bpnreg` is the only
R-package providing Bayesian regression and mixed-effects models for
circular outcomes based on the projected normal distribution.

From its first release in 2018 at least eight published articles
[@Tyson-Carr2020-mu; @Cremers2021-mm; @Cote2020-xg; @Rafferty2020-pb; @Ojeda_undated-by; @Olson2020-al; @Klugkist2018-ag; @Spinks2019-ya]
have used `bpnreg` in their analyses.

# Example Usage

To use `bpnreg` the user first needs install the package from `CRAN` and
load the package as follows:

``` {.r}
install.packages("bpnreg")
library(bpnreg)
```

    ## Warning: package 'bpnreg' was built under R version 4.0.3

The two main functions are `bpnr` for fitting Bayesian projected normal
multiple regression models for circular outcomes and `bpnme` for fitting
Bayesian projected normal mixed-effects models for circular outcomes.

For more examples we refer to @Cremers2018-ys. Answers to frequently
asked questions regarding the use of `bpnreg` can be found in the `FAQ`
vignette available on `CRAN`.

# Acknowledgements

JC is supported for this work by a research grant from the Novo Nordisk
Foundation ("Harnessing The Power of Big Data to Address the Societal
Challenge of Aging." NNF17OC0027812)

# References {#references .unnumbered}