#' bpnreg: A package to analyze Bayesian projected normal circular regression models
#'
#' This package contains functions to analyze circular regression type models
#' (multivariate and mixed-effects). It is based on the 'embedding' approach to
#' circular modelling and makes use of the projected normal distribution. Its
#' estimation method is a Bayesian MCMC sampler. Further technical details can
#' be found in Cremers, Mulder & Klugkist (2017) and Cremers & Klugkist (2017,
#' under review).
#'
#' A tutorial on how to use this package can be found in Cremers & Klugkist
#' (2017, working paper). More details on the sampling algorithm and
#' interpretation of the coefficients from the model can be found in Cremers,
#' Mulder & Klugkist (2017) and Cremers, Mainhard & Klugkist (2017, under review).
#'
#' @section Functions: The main functions of the package are:
#'
#'   \code{\link{bpnr}}, which runs an MCMC sampler in \code{C++} and
#'   returns an S3 object of type \code{bpnr}, which can be further
#'   analyzed through associated functions.
#'
#'   and
#'
#'   \code{\link{bpnme}}, which runs an MCMC sampler in \code{R} and
#'   returns an S3 object of type \code{bpnme}, which can be further
#'   analyzed through associated functions.
#'
#' @section Datasets: Datasets included in this package are:
#'
#'   \code{\link{Motor}}, A dataset from a research by Puglisi et.al. (2017) on
#'   the role of attention in human motor resonance.
#'
#'   and
#'
#'   \code{\link{Maps}}, A dataset from a research by Warren et.al. (2017) on
#'   the geometry of humans' knowledge of navigation space.
#'
#' @source Cremers, J., Mulder, K.T. & Klugkist, I. (In press). Circular
#'   interpretation of regression coefficients. British Journal of Mathematical
#'   and Statistical Psychology.
#'
#' @source Cremers, J., Mainhard, M.T. & Klugkist, I. (2017). Assessing a Bayesian
#'   Embedding Approach to Circular Regression Models. Manuscript under review.
#'
#' @source Cremers, J., & Klugkist, I. (2017). How to analyze circular data: A
#'    tutorial for projected normal regression models. Working paper.
#'
#' @useDynLib bpnreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd rnorm runif pnorm dnorm var model.matrix predict rWishart terms
#' @importFrom stats model.frame plot.ts reformulate as.formula formula
#' @importFrom MASS mvrnorm
#' @importFrom haven read_spss
#' @importFrom methods is
#' @importFrom utils combn
#'
#' @docType package
#' @name bpnreg

NULL
