#' Fit a Bayesian circular regression model
#'
#' This function fits a Bayesian circular regression model based on the
#' projected normal distribution.
#'
#' @param pred.I model equation for effects of component 1.
#' @param data the dataframe used for analysis.
#' @param pred.II model equation for effects of component 2.
#' @param its output iterations of the MCMC sampler.
#' @param burn number of burn-in iterations.
#' @param n.lag amount of lag for the iterations and burn-in.
#' @param seed user-specified random seed.
#'
#' @details Because the model is based on the projected normal distribution, a
#'   model equation has to be given for two components. By default the equation
#'   of the second component \code{pred.II} is set to be equal to that of the
#'   first component. Note that the circular outcome needs to be measured in
#'   radians on a scale from 0 to 2\eqn{\pi}. For more information about the
#'   projected normal distribution see Presnell, Morrisson & Littell (1998).
#'
#'   A tutorial on how to use this function can be found in Cremers & Klugkist
#'   (2018). More details on the sampling algorithm and interpretation of the
#'   coefficients from the model can be found in Cremers, Mulder & Klugkist
#'   (2018) and Cremers, Mainhard & Klugkist (2018). The uninformative priors
#'   for the regression coefficients of the two components are set to N(0, 10000).
#'
#' @return A \code{bpnr} object, which can be further analyzed using the
#'   associated functions \code{\link{traceplot.bpnr}}, \code{\link{BFc.bpnr}},
#'   \code{\link{coef_lin.bpnr}}, \code{\link{coef_circ.bpnr}},
#'   \code{\link{fit.bpnr}} and \code{\link{print.bpnr}}.
#'
#'   A \code{bpnr} object contains the following elements (some elements are not
#'   returned if not applicable)
#'
#'   \describe{ \item{\code{beta1}}{A matrix of posterior samples for the
#'   coefficients \code{beta1} of the first component.} \item{\code{beta2}}{A matrix
#'   of posterior samples for the coefficients \code{beta2} for the second
#'   component.} \item{\code{Likelihood}}{A matrix containing the posterior
#'   density values for all individuals in the dataset for all iterations. The
#'   rowsums of this matrix are the likelihood values for all iterations}
#'   \item{\code{its}}{Number of output iterations.} \item{\code{n.lag}}{One in
#'   \code{n.lag} iterations will be saved as output iteration. Set lag to 1 to
#'   save all iterations (default).} \item{\code{burn-in}}{Burn-in time for the
#'   MCMC sampler.} \item{\code{p1}}{Number of parameters predicting the first
#'   component.} \item{\code{p2}}{Number of parameters predicting the second
#'   component.} \item{\code{theta}}{The circular outcome vector measured in radians.}
#'   \item{\code{a.x}}{A matrix of posterior samples for \code{a.x} which
#'   describes the location of the inflection point of the regression curve on
#'   the axis of the predictor.} \item{\code{a.c}}{A matrix of posterior samples
#'   for \code{a.c} which describes the location of the inflection point of the
#'   regression curve on the axis of the circular outcome.} \item{\code{b.c}}{A
#'   matrix of posterior samples for \code{b.c} which describes the slope of the
#'   tangent line at the inflection point.} \item{\code{SAM}}{A matrix of
#'   posterior samples for the circular regression slopes at the mean.}
#'   \item{\code{AS}}{A matrix of posterior samples for the average slopes of
#'   the circular regression.} \item{\code{SSDO}}{A matrix of posterior samples
#'   for the signed shortest distance to the origin.} \item{\code{circ.diff}}{A
#'   matrix of posterior samples for the circular difference between levels of
#'   categorical variables and the intercept.} \item{\code{Call}}{The matched
#'   call.} \item{\code{lin.coef.I}}{The mean, mode, standard deviation and 95 %
#'   confidence interval of the highest posterior density of the linear
#'   coefficients for \code{beta1}.} \item{\code{lin.coef.II}}{The mean, mode,
#'   standard deviation and 95 % confidence interval of the highest posterior
#'   density of the linear coefficients for \code{beta2}.}
#'   \item{\code{circ.coef}}{The mean, mode, standard deviation and 95 %
#'   confidence interval of the highest posterior density for the \code{a.x},
#'   \code{a.c}, \code{b.c}, \code{AS}, \code{SAM} and \code{SSDO} of the
#'   circular coefficients.} \item{\code{circ.coef.cat}}{The mean, mode,
#'   standard deviation and 95 % confidence interval of the highest posterior
#'   density the circular difference between levels of categorical variables and
#'   the intercept.} \item{\code{circ.coef.means}}{The mean, mode, standard
#'   deviation and 95 % confidence interval of the highest posterior density of
#'   circular means of the categorical variables.} \item{\code{model.fit}}{A
#'   list of information criteria for assessment of model fit.}
#'   \item{\code{mm}}{A list of information, model matrices, sample size, etc.
#'   on the specified model.} }
#'
#'
#' @source Cremers, J., Mulder, K.T. & Klugkist, I. (2018). Circular
#'   interpretation of regression coefficients. British Journal of Mathematical
#'   and Statistical Psychology, 71(1), 75-95.
#'
#' @source Cremers, J., Mainhard, M.T. & Klugkist, I. (2018). Assessing a
#'   Bayesian Embedding Approach to Circular Regression Models. Methodology, 14,
#'   69-81.
#'
#' @source Cremers, J. & Klugkist, I. (2018). One direction? A tutorial for
#'   circular data with examples in cognitive psychology. Frontiers in
#'   Psychology: Cognitive Science.
#'
#' @source Presnell, B., Morrison, S.P. & Littell, R.C. (1998). Projected
#'   multivariate linear models for directional data. Journal of the American
#'   Statistical Association, 93 (443), 1068 - 1077.
#'
#' @examples
#' library(bpnreg)
#' bpnr(Phaserad ~ Cond + AvAmp, Motor)
#'
#' @export
#'

bpnr <- function(pred.I, data, pred.II = pred.I,
                 its = 1000, burn = 1, n.lag = 1,
                 seed = NULL){

  if (!is.null(seed)){set.seed(seed)}

  mm <- mmr(pred.I, data, pred.II)

  output <- pnr(mm$theta, mm$XI, mm$XII, its, n.lag, burn)

  class(output) <- c("bpnr", class(output))

  summary.stats <- sumr(output, mm)

  output$beta1 <- summary.stats$beta1
  output$beta2 <- summary.stats$beta2
  output$a.x <- summary.stats$a.x
  output$a.c <- summary.stats$a.c
  output$b.c <- summary.stats$b.c
  output$SAM <- summary.stats$SAM
  output$AS <- summary.stats$AS
  output$SSDO <- summary.stats$SSDO
  output$circ.diff <- summary.stats$circ.diff

  output$Call <- match.call()
  output$lin.coef.I <- summary.stats$lin.res.I
  output$lin.coef.II <- summary.stats$lin.res.II
  output$circ.coef <- summary.stats$circ.res
  output$circ.coef.cat <- summary.stats$circ.res.cat
  output$circ.coef.means <- summary.stats$circ.res.means
  output$model.fit <- summary.stats$model.fit
  output$mm <- mm

  output

}

#' Fit a Bayesian circular mixed-effects model
#'
#' This function fits a Bayesian circular mixed-effects model based on the
#' projected normal distribution.
#'
#' @param pred.I model equation for effects of component 1.
#' @param data the dataframe used for analysis.
#' @param pred.II model equation for effects of component 2.
#' @param its output iterations of the MCMC sampler.
#' @param burn number of burn-in iterations.
#' @param n.lag amount of lag for the iterations and burn-in.
#' @param seed user-specified random seed.
#'
#' @details Because the model is based on the projected normal distribution, a
#'   model equation has to be given for the fixed and random effects of the two
#'   components. By default the model equation of the second component
#'   \code{pred.II} is set to be equal to that of the first component. Note that
#'   the circular outcome needs to be measured in radians on a scale from 0 to
#'   2\eqn{\pi}. For more information about the projected normal distribution
#'   see Presnell, Morrisson & Littell (1998). The model can handle at most one
#'   grouping factor.
#'
#'   A tutorial on how to use this function can be found in Cremers & Klugkist
#'   (2018). More details on the sampling algorithm and interpretation of the
#'   coefficients from the model can be found in Nuñez-Antonio & Guttiérrez-Peña
#'   (2014) and Cremers, Pennings, Mainhard & Klugkist (2019). The uninformative priors
#'   for the fixed effect regression coefficients of the two components are set to N(0, 10000).
#'   Note that the model is only developed for models with a single nesting variable.
#'
#' @return A \code{bpnme} object, which can be further analyzed using the
#'   associated functions \code{\link{traceplot.bpnme}},
#'   \code{\link{BFc.bpnme}}, \code{\link{coef_lin.bpnme}},
#'   \code{\link{coef_circ.bpnme}}, \code{\link{coef_ran.bpnme}},
#'   \code{\link{fit.bpnme}} and \code{\link{print.bpnme}}.
#'
#'   A \code{bpnr} object contains the following elements (some elements are not
#'   returned if not applicable)
#'
#'   \describe{ \item{\code{beta1}}{A matrix of posterior samples for the fixed
#'   effects coefficients for the first component.} \item{\code{beta2}}{A
#'   matrix of posterior samples for the fixed effects coefficients for the
#'   second component.} \item{\code{b1}}{An array of posterior samples for the
#'   random effects coefficients for the first component.} \item{\code{b2}}{An
#'   array of posterior samples for the random effects coefficients for the
#'   second component.} \item{\code{omega1}}{An array of posterior samples for
#'   the random effect variances of the first component.}
#'   \item{\code{omega2}}{An array of posterior samples for the random effect
#'   variances of the second component.} \item{\code{predictiva}}{A list
#'   containing the posterior density values for all timepoints of individuals
#'   in the dataset for all iterations. The rowsums of this matrix are the
#'   likelihood values for all iterations} \item{\code{circular.ri}}{A vector of
#'   posterior samples for the circular random intercepts.}
#'   \item{\code{N}}{Number of observed cases.} \item{\code{its}}{Number of
#'   output iterations.} \item{\code{n.lag}}{One in \code{n.lag} iterations will
#'   be saved as output iteration. Set lag to 1 to save all iterations
#'   (default).} \item{\code{burn}}{Burn-in time for the MCMC sampler.}
#'   \item{\code{p1}}{Number of fixed effect parameters predicting the first
#'   component.} \item{\code{p2}}{Number of fixed effect parameters predicting
#'   the second component.} \item{\code{q1}}{Number of random effect parameters
#'   predicting the first component.} \item{\code{q2}}{Number of random effect
#'   parameters predicting the second component.} \item{\code{a.x}}{A matrix of
#'   posterior samples for \code{a.x} which describes the location of the
#'   inflection point of the regression curve on the axis of the predictor.}
#'   \item{\code{a.c}}{A matrix of posterior samples for \code{a.c} which
#'   describes the location of the inflection point of the regression curve on
#'   the axis of the circular outcome.} \item{\code{b.c}}{A matrix of posterior
#'   samples for \code{b.c} which describes the slope of the tangent line at the
#'   inflection point.} \item{\code{SAM}}{A matrix of posterior samples for the
#'   circular regression slopes at the mean.} \item{\code{AS}}{A matrix of
#'   posterior samples for the average slopes of the circular regression.}
#'   \item{\code{SSDO}}{A matrix of posterior samples for the signed shortest
#'   distance to the origin.} \item{\code{circ.diff}}{A matrix of posterior
#'   samples for the circular difference found between levels of categorical
#'   variables and the intercept.} \item{\code{cRSnum}}{A string indicating
#'   whether there are continuous variables with a random slope}
#'   \item{\code{cRScat}}{A string indicating whether there are categorical
#'   variables with a random slope} \item{\code{cRS}}{A string indicating
#'   whether there are categorical or continuous variables with a random slope}
#'   \item{\code{cRI}}{A vector of posterior samples of the mean resultant
#'   length of the circular random intercept, a measure of concentration.}
#'   \item{\code{Call}}{The matched call.} \item{\code{lin.coef.I}}{The mean,
#'   mode, standard deviation and 95 % confidence interval of the highest
#'   posterior density of the linear fixed effect coefficients for \code{beta1}.}
#'   \item{\code{lin.coef.II}}{The mean, mode, standard deviation and 95 %
#'   confidence interval of the highest posterior density of the linear fixed
#'   effect coefficients for \code{beta2}.} \item{\code{circ.coef}}{The mean, mode,
#'   standard deviation and 95 % confidence interval of the highest posterior
#'   density for \code{a.x}, \code{a.c}, \code{SSDO}, and the circular fixed
#'   effect coefficients \code{b.c}, \code{AS}, and \code{SAM}}
#'   \item{\code{circ.coef.cat}}{The mean, mode, standard deviation and 95 %
#'   confidence interval of the highest posterior density the circular
#'   difference between levels of categorical variables and the intercept.}
#'   \item{\code{circ.coef.means}}{The mean, mode, standard deviation and 95 %
#'   confidence interval of the highest posterior density of circular means of
#'   the categorical variables.} \item{\code{model.fit}}{A list of information
#'   criteria for assessment of model fit.} \item{\code{lin.res.varrand.I}}{The
#'   mean, mode, standard deviation and 95 % confidence interval of  the
#'   variances of the random intercepts and slopes of component I.}
#'   \item{\code{lin.res.varrand.II}}{The mean, mode, standard deviation and 95
#'   % confidence interval of the variances of the random intercepts and slopes
#'   of component II.} \item{\code{circ.res.varrand}}{The mean, mode, standard
#'   deviation and 95 % confidence interval of the circular variances of the
#'   random intercepts and slopes.} \item{\code{mm}}{A list of information,
#'   model matrices, sample size, etc. on the specified model.} }
#'
#' @source Cremers, J., Mainhard, M.T. & Klugkist, I. (2018). Assessing a
#'   Bayesian Embedding Approach to Circular Regression Models. Methodology, 14, 69-81.
#'
#' @source Cremers, J. & Klugkist, I. (2018). One direction? A tutorial for
#'   circular data with examples in cognitive psychology. Frontiers in
#'   Psychology: Cognitive Science.
#'
#' @source Cremers, J., Pennings, H.J.M., Mainhard, M.T. & Klugkist, I. (2019).
#'   Circular Modelling of Circumplex Measurements for Interpersonal Behavior.
#'   Assessment, Online First.
#'
#' @source Nuñez-Antonio, G. & Gutiérrez-Peña, E. (2014). A Bayesian model for
#'   longitudinal circular data based on the projected normal distribution.
#'   Computational Statistics and Data Analysis, 71, 506-519.
#'
#' @source Presnell, B., Morrison, S.P. & Littell, R.C. (1998). Projected
#'   multivariate linear models for directional data. Journal of the American
#'   Statistical Association, 93 (443), 1068 - 1077.
#'
#' @examples
#' library(bpnreg)
#' bpnme(Error.rad ~ Maze + Trial.type + (1|Subject), Maps, its = 100)
#'
#' @export

bpnme <- function(pred.I, data, pred.II = pred.I,
                  its = 1000, burn = 1, n.lag = 1,
                  seed = NULL){

  if (!is.null(seed)){set.seed(seed)}

  mm <- mmme(pred.I, data, pred.II)

  circular.ri <- matrix(NA, mm$N, its)

  predictiva <- list()

  for (i in 1:mm$N){
    predictiva[[i]] <- matrix(NA, its, mm$no.Meas[i])
  }


  output <- pnme(lapply(mm$theta, cos), lapply(mm$theta, sin),
                 mm$XI, mm$XII, mm$ZI, mm$ZII, mm$ZtZI, mm$ZtZII,
                 mm$R, predictiva,
                 its = its, burn = burn, lag = n.lag, mm$N)

  class(output) <- c("bpnme", class(output))

  rownames(output$beta1) <- colnames(mm$XI[[1]])
  rownames(output$beta2) <- colnames(mm$XII[[1]])
  colnames(output$b1) <- colnames(mm$ZI[[1]])
  colnames(output$b2) <- colnames(mm$ZII[[1]])
  colnames(output$omega1) <- colnames(mm$ZI[[1]])
  rownames(output$omega1) <- colnames(mm$ZI[[1]])
  colnames(output$omega2) <- colnames(mm$ZII[[1]])
  rownames(output$omega2) <- colnames(mm$ZII[[1]])

  # Compute circular RI

  for (i in 1:mm$N){

    if(!("(Intercept)" %in% colnames(mm$XI[[i]])) & !("(Intercept)" %in% colnames(mm$XII[[i]]))){

      circular.ri[i,] <- atan2(0 + output$b2[i,1,], 0 + output$b1[i,1,])

    }else if(!("(Intercept)" %in% colnames(mm$XI[[i]]))){

      circular.ri[i,] <- atan2(output$beta2[1,] + output$b2[i,1,], 0 + output$b1[i,1,])

    }else if(!("(Intercept)" %in% colnames(mm$XII[[i]]))){

      circular.ri[i,] <- atan2(0 + output$b2[i,1,], output$beta1[1,] + output$b1[i,1,])

    }else{

      circular.ri[i,] <- atan2(output$beta2[1,] + output$b2[i,1,], output$beta1[1,] + output$b1[i,1,])

    }
  }

  output[["circular.ri"]] = circular.ri
  output$beta1 = t(output$beta1)
  output$beta2 = t(output$beta2)

  summary.stats <- summe(output, mm)

  output$a.x <- summary.stats$a.x
  output$a.c <- summary.stats$a.c
  output$b.c <- summary.stats$b.c
  output$SAM <- summary.stats$SAM
  output$AS <- summary.stats$AS
  output$SSDO <- summary.stats$SSDO
  output$circ.diff <- summary.stats$circ.diff
  output$cRSnum <- summary.stats$varrand.num
  output$cRScat <- summary.stats$varrand.cat
  output$cRS <- cbind(summary.stats$varrand.cat, summary.stats$varrand.num)
  output$cRI <- summary.stats$circ.varrand.ri

  output$Call <- match.call()
  output$lin.coef.I <- summary.stats$lin.res.I
  output$lin.coef.II <- summary.stats$lin.res.II
  output$circ.coef <- summary.stats$circ.res
  output$circ.coef.cat <- summary.stats$circ.res.cat
  output$circ.coef.means <- summary.stats$circ.res.means
  output$model.fit <- summary.stats$model.fit
  output$lin.res.varrand.I <- summary.stats$lin.res.varrand.I
  output$lin.res.varrand.II <- summary.stats$lin.res.varrand.II
  output$circ.res.varrand <- summary.stats$circ.res.varrand
  output$mm <- mm

  output

}

