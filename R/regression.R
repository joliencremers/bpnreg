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
#'   first component. For more information about the projected normal
#'   distribution see Presnell, Morrisson & Littell (1998).
#'
#'   A tutorial on how to use this function can be found in Cremers & Klugkist
#'   (2017, working paper). More details on the sampling algorithm and
#'   interpretation of the coefficients from the model can be found in Cremers,
#'   Mulder & Klugkist (2018) and Cremers, Mainhard & Klugkist (in press).
#'
#' @return A \code{bpnr} object, which can be further analyzed using the
#'   associated functions \code{\link{traceplot.bpnr}}, \code{\link{BFc.bpnr}},
#'   \code{\link{coef_lin.bpnr}}, \code{\link{coef_circ.bpnr}},
#'   \code{\link{residuals.bpnr}} ,\code{\link{predict.bpnr}},
#'   \code{\link{fit.bpnr}} and \code{\link{print.bpnr}}.
#'
#'   A \code{bpnr} object contains the following elements (some elements are not
#'   returned if not applicable)
#'
#'   \describe{
#'   \item{\code{B1}}{A matrix of posterior samples for the coefficients \code{B1}
#'   of the first component.}
#'   \item{\code{B2}}{A matrix of posterior samples for the coefficients \code{B2}
#'   for the second component.}
#'   \item{\code{Likelihood}}{A matrix containing the posterior density values for all individuals in the dataset
#'    for all iterations. The rowsums of this matrix are the likelihood values for all iterations}
#'   \item{\code{its}}{Number of output iterations.}
#'   \item{\code{n.lag}}{One in \code{n.lag} iterations will be saved as output iteration. Set lag to 1 to save all iterations (default).}
#'   \item{\code{burn-in}}{Burn-in time for the MCMC sampler.}
#'   \item{\code{p1}}{Number of parameters predicting the first component.}
#'   \item{\code{p2}}{Number of parameters predicting the second component.}
#'   \item{\code{theta}}{The circular outcome vector.}
#'   \item{\code{a.x}}{A matrix of posterior samples for \code{a.x} which describes the location of the
#'   inflection point of the regression curve on the axis of the predictor.}
#'   \item{\code{a.c}}{A matrix of posterior samples for \code{a.c} which describes the location of the
#'   inflection point of the regression curve on the axis of the circular outcome.}
#'   \item{\code{b.c}}{A matrix of posterior samples for \code{b.c} which describes the slope of the tangent line at the inflection point.}
#'   \item{\code{SAM}}{A matrix of posterior samples for the circular regression slopes at the mean.}
#'   \item{\code{AS}}{A matrix of posterior samples for the average slopes of the circular regression.}
#'   \item{\code{SSDO}}{A matrix of posterior samples for the signed shortest distance to the origin.}
#'   \item{\code{circ.diff}}{A matrix of posterior samples for the circular difference between levels of categorical variables and the intercept.}
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{lin.coef.I}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of the linear coefficients for \code{B1}.}
#'   \item{\code{lin.coef.II}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of the linear coefficients for \code{B2}.}
#'   \item{\code{circ.coef}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density for the \code{a.x}, \code{a.c},
#'   \code{b.c}, \code{AS}, \code{SAM} and \code{SSDO} of the circular coefficients.}
#'   \item{\code{circ.coef.cat}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density the circular difference between levels of categorical variables and the intercept.}
#'   \item{\code{circ.coef.means}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of circular means of the categorical variables.}
#'   \item{\code{model.fit}}{A list of information criteria for assesment of model fit.}
#'   \item{\code{mm}}{A list of information, model matrices, sample size, etc. on the specified model.}
#'   }
#'
#'
#' @source Cremers, J., Mulder, K.T. & Klugkist, I. (2018). Circular
#'   interpretation of regression coefficients. British Journal of Mathematical
#'   and Statistical Psychology, 71(1), 75-95.
#'
#' @source Cremers, J., Mainhard, M.T. & Klugkist, I. (in press). Assessing a
#'   Bayesian Embedding Approach to Circular Regression Models. Methodology.
#'
#' @source Cremers, J., & Klugkist, I. (2017). How to analyze circular data: A
#'   tutorial for projected normal regression models. Under review.
#'
#' @source Presnell, B., Morrison, S.P. & Littell, R.C. (1998). Projected
#'   multivariate linear models for directional data. Journal of the Americal
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

  output$B1 <- summary.stats$B1
  output$B2 <- summary.stats$B2
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

#' Fit a Baysian circular mixed-effects model
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
#'   model equation has to be given for the fixed and random effects of thetwo
#'   components. By default the model equation of the second component
#'   \code{pred.II} is set to be equal to that of the
#'   first component. For more information about the projected normal
#'   distribution see Presnell, Morrisson & Littell (1998).
#'
#'   A tutorial on how to use this function can be found in Cremers & Klugkist
#'   (2017, working paper). More details on the sampling algorithm and
#'   interpretation of the coefficients from the model can be found in Cremers & Klugkist (2017, working paper).
#'
#' @return A \code{bpnme} object, which can be further analyzed using the
#'   associated functions \code{\link{traceplot.bpnme}}, \code{\link{BFc.bpnme}},
#'   \code{\link{coef_lin.bpnme}}, \code{\link{coef_circ.bpnme}}, \code{\link{coef_ran.bpnme}},
#'   \code{\link{residuals.bpnme}}, \code{\link{predict.bpnme}},
#'   \code{\link{fit.bpnme}} and \code{\link{print.bpnme}}.
#'
#'   A \code{bpnr} object contains the following elements (some elements are not
#'   returned if not applicable)
#'
#'   \describe{
#'   \item{\code{Beta.I}}{A matrix of posterior samples for the fixed effects coefficients for the first component.}
#'   \item{\code{Beta.II}}{A matrix of posterior samples for the fixed effects coefficients for the second component.}
#'   \item{\code{B.I}}{An array of posterior samples for the random effects coefficients for the first component.}
#'   \item{\code{B.II}}{An array of posterior samples for the random effects coefficients for the second component.}
#'   \item{\code{VCovI}}{An array of posterior samples for the random effect variances of the first component.}
#'   \item{\code{VCovII}}{An array of posterior samples for the random effect variances of the second component.}
#'   \item{\code{predictiva}}{A list containing the posterior density values for all timepoints of individuals in the dataset
#'    for all iterations. The rowsums of this matrix are the likelihood values for all iterations}
#'   \item{\code{circular.ri}}{A vector of posterior samples for the circular random intercepts.}
#'   \item{\code{N}}{Number of observed cases.}
#'   \item{\code{its}}{Number of output iterations.}
#'   \item{\code{n.lag}}{One in \code{n.lag} iterations will be saved as output iteration. Set lag to 1 to save all iterations (default).}
#'   \item{\code{burn}}{Burn-in time for the MCMC sampler.}
#'   \item{\code{p1}}{Number of fixed effect parameters predicting the first component.}
#'   \item{\code{p2}}{Number of fixed effect parameters predicting the second component.}
#'   \item{\code{q1}}{Number of random effect parameters predicting the first component.}
#'   \item{\code{q2}}{Number of random effect parameters predicting the second component.}
#'   \item{\code{a.x}}{A matrix of posterior samples for \code{a.x} which describes the location of the
#'   inflection point of the regression curve on the axis of the predictor.}
#'   \item{\code{a.c}}{A matrix of posterior samples for \code{a.c} which describes the location of the
#'   inflection point of the regression curve on the axis of the circular outcome.}
#'   \item{\code{b.c}}{A matrix of posterior samples for \code{b.c} which describes the slope of the tangent line at the inflection point.}
#'   \item{\code{SAM}}{A matrix of posterior samples for the circular regression slopes at the mean.}
#'   \item{\code{AS}}{A matrix of posterior samples for the average slopes of the circular regression.}
#'   \item{\code{SSDO}}{A matrix of posterior samples for the signed shortest distance to the origin.}
#'   \item{\code{circ.diff}}{A matrix of posterior samples for the circular difference found between levels of categorical variables and the intercept.}
#'   \item{\code{cRSnum}}{A string indicating whether there are continuous variables with a random slope}
#'   \item{\code{cRScat}}{A string indicating whether there are categorical variables with a random slope}
#'   \item{\code{cRS}}{A string indicating whether there are categorical or continuous variables with a random slope}
#'   \item{\code{cRI}}{A vector of posterior samples of the mean resultant length of the circular random intercept, a measure of concentration.}
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{lin.coef.I}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of the linear fixed effect coefficients for \code{B1}.}
#'   \item{\code{lin.coef.II}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of the linear fixed effect coefficients for \code{B2}.}
#'   \item{\code{circ.coef}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density for \code{a.x}, \code{a.c}, \code{SSDO}, and the circular fixed effect coefficients
#'   \code{b.c}, \code{AS}, and \code{SAM}}
#'   \item{\code{circ.coef.cat}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density the circular difference between levels of categorical variables and the intercept.}
#'   \item{\code{circ.coef.means}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the highest posterior density of circular means of the categorical variables.}
#'   \item{\code{model.fit}}{A list of information criteria for assesment of model fit.}
#'   \item{\code{lin.res.varrand.I}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of  the variances of the random intercepts and slopes of component I.}
#'   \item{\code{lin.res.varrand.II}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the variances of the random intercepts and slopes of component II.}
#'   \item{\code{circ.res.varrand}}{The mean, mode, standard deviation and 95 % confidence
#'   interval of the circular variances of the random intercepts and slopes.}
#'   \item{\code{mm}}{A list of information, model matrices, sample size, etc. on the specified model.}
#'   }
#'
#' @source Cremers, J., Mainhard, M.T. & Klugkist, I. (in press). Assessing a
#'   Bayesian Embedding Approach to Circular Regression Models. Methodology
#'
#' @source Cremers, J., & Klugkist, I. (2017). How to analyze circular data: A
#'   tutorial for projected normal regression models. Under review.
#'
#' @source Cremers, J., & Klugkist, I. (2017). Longitudinal
#'   circular modelling of circumplex measurements for teacher behavior. Working
#'   paper.
#'
#' @source Presnell, B., Morrison, S.P. & Littell, R.C. (1998). Projected
#'   multivariate linear models for directional data. Journal of the Americal
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

  if(!"(Intercept)" %in% colnames(mm$mm_ran.I) | !"(Intercept)" %in% colnames(mm$mm_ran.II)){

    stop("No random intercept in the model")

  }

  if(!all(colnames(mm$mm_ran.I) %in% colnames(mm$mm.I)) | !all(colnames(mm$mm_ran.II) %in% colnames(mm$mm.II))){

    stop("Not all random effects have a corresponding fixed effect")

  }

  burn <- burn * n.lag
  tm <- burn + (its * n.lag)
  p1 <- length(mm$XI[[mm$N]][1, ])
  p2 <- length(mm$XII[[mm$N]][1, ])
  q1 <- length(mm$ZI[[mm$N]][1, ])
  q2 <- length(mm$ZII[[mm$N]][1, ])
  OmegaI  <- diag(q1)
  OmegaII <- diag(q2)

  #priors

  v1  <- q1
  v2  <- q2

  B1  <- diag(v1) * 0.001
  B2  <- diag(v2) * 0.001

  A1  <- matrix(0, p1, p1)
  A2  <- matrix(0, p2, p2)

  #matrices for final results
  Beta.I   <- matrix(NA, its, p1)
  Beta.II  <- matrix(NA, its, p2)
  colnames(Beta.I) <- colnames(mm$XI[[1]])
  colnames(Beta.II) <- colnames(mm$XII[[1]])
  B.I      <- array(NA, c(mm$N, q1, its))
  B.II     <- array(NA, c(mm$N, q2, its))
  colnames(B.I) <- colnames(mm$ZI[[1]])
  colnames(B.II) <- colnames(mm$ZII[[1]])
  VCovI    <- array(NA, dim = c(q1, q1, its))
  VCovII   <- array(NA, dim = c(q2, q2, its))
  colnames(VCovI) <- colnames(mm$ZI[[1]])
  colnames(VCovII) <- colnames(mm$ZII[[1]])
  rownames(VCovI) <- colnames(mm$ZI[[1]])
  rownames(VCovII) <- colnames(mm$ZII[[1]])
  int_rand_c <- c()
  circular.ri <- matrix(NA, mm$N, its)


  Y1       <- list()
  Y2       <- list()
  YI       <- list()
  YII      <- list()
  predictiva <- list()

  for (i in 1:mm$N){
    predictiva[[i]] <- matrix(NA, its, mm$no.Meas[i])
  }


  if (q1 == 1 & q2 == 1){

    ####real iterations####
    for (k in 1:(tm)){

      for (i in 1:mm$N){
        YI[[i]]  <- mm$R[[i]] * cos(mm$theta[[i]])
        YII[[i]] <- mm$R[[i]] * sin(mm$theta[[i]])
      }

      #simulations for beta, bi's and Omega

      beta.I    <- betaBlock.fRI(OmegaI, YI, mm$XI, mm$ZI, p1, A1, mm$N)
      beta.II   <- betaBlock.fRI(OmegaII, YII, mm$XII, mm$ZII, p2, A2, mm$N)
      b.I       <- b.f(OmegaI, beta.I, YI, mm$XI, q1, mm$ZI, mm$ZtZI, mm$N)
      b.II      <- b.f(OmegaII, beta.II, YII, mm$XII, q2, mm$ZII, mm$ZtZII, mm$N)
      OmegaI    <- Omega.f(b.I, B1, v1, q1, mm$N)
      OmegaII   <- Omega.f(b.II, B2, v2, q2, mm$N)

      for (i in 1:mm$N){

        if(!("(Intercept)" %in% colnames(mm$XI[[i]])) & !("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i], 0+b.I[i])

        }else if(!("(Intercept)" %in% colnames(mm$XI[[i]]))){

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i,1], 0+b.I[i])

        }else if(!("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i], beta.I[1]+b.I[i])

        }else{

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i], beta.I[1]+b.I[i])

        }

        for (j in 1:mm$no.Meas[i]){

          t.aux     <- mm$theta[[i]][j]
          mu.ij.I   <- c(beta.I %*% mm$XI[[i]][j, ] + b.I[i])
          mu.ij.II  <- c(beta.II %*% mm$XII[[i]][j, ] + b.II[i])
          bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
          mm$R[[i]][j] <- slice_r_me(t.aux, mu.ij.I, mu.ij.II, mm$R[[i]][j])

          if (k - burn > 0){
            if ( ( (k - burn) %% n.lag) == 0){
              ii <- (k - burn) / n.lag

              norm2  <- mu.ij.I^2 + mu.ij.II^2
              c         <- (1 + ( (bb * pnorm(bb)) / dnorm(bb)))

              predictiva[[i]][ii,j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2))*c
            }
          }
        }
      }

      if (k - burn <= 0 & ( (burn-k) / n.lag) %% 50 == 0){

        ii <- (burn - k) / n.lag; print(paste("burn-in iteration", ii))

        }

      #save values for each iteration
      if(k - burn > 0){
        if ( ( (k-burn) %% n.lag) == 0){

          ii <- (k - burn) / n.lag
          if ( (ii %% 50) == 0) {print(paste("iteration", ii))}

          Beta.I[ii, ]  <- beta.I
          Beta.II[ii, ] <- beta.II
          B.I[,,ii]    <- b.I
          B.II[,,ii]   <- b.II
          VCovI[,,ii]  <- solve(OmegaI)
          VCovII[,,ii] <- solve(OmegaII)
          Y1[[ii]] <- YI
          Y2[[ii]] <- YII
          circular.ri[, ii] <- int_rand_c

        }
      }
    }




  }else if(q1 >= 1 & q2 == 1){


    ####real iterations####
    for(k in 1:(tm))
    {

      for(i in 1:mm$N){
        YI[[i]]  <- mm$R[[i]] * cos(mm$theta[[i]])
        YII[[i]] <- mm$R[[i]] * sin(mm$theta[[i]])
      }

      #simulations for beta, bi's and Omega
      beta.I    <- betaBlock.fRS(OmegaI, YI, mm$XI, mm$ZI, p1, A1, mm$N)
      beta.II   <- betaBlock.fRI(OmegaII, YII, mm$XII, mm$ZII, p2, A2, mm$N)
      b.I       <- b.f(OmegaI, beta.I, YI, mm$XI, q1, mm$ZI, mm$ZtZI, mm$N)
      b.II      <- b.f(OmegaII, beta.II, YII, mm$XII, q2, mm$ZII, mm$ZtZII, mm$N)
      OmegaI    <- Omega.f(b.I, B1, v1, q1, mm$N)
      OmegaII   <- Omega.f(b.II, B2, v2, q2, mm$N)

      for(i in 1:mm$N){

        if(!("(Intercept)" %in% colnames(mm$XI[[i]])) & !("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i], 0+b.I[i,1])

        }else if(!("(Intercept)" %in% colnames(mm$XI[[i]]))){

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i], 0+b.I[i,1])

        }else if(!("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i], beta.I[1]+b.I[i,1])

        }else{

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i], beta.I[1]+b.I[i,1])

        }

        for(j in 1:mm$no.Meas[i]){

          t.aux     <- mm$theta[[i]][j]
          mu.ij.I   <- c(beta.I %*% mm$XI[[i]][j, ]+b.I[i, ]%*%mm$ZI[[i]][j, ])
          mu.ij.II  <- c(beta.II %*% mm$XII[[i]][j, ] + b.II[i])
          bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
          mm$R[[i]][j] <- slice_r_me(t.aux, mu.ij.I, mu.ij.II, mm$R[[i]][j])

          if (k - burn > 0){
            if ( ( (k - burn) %% n.lag) == 0){
              ii <- ( k - burn) / n.lag

              norm2  <- mu.ij.I^2 + mu.ij.II^2
              c         <- (1 + ( (bb * pnorm(bb)) / dnorm(bb)))

              predictiva[[i]][ii,j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2))*c
            }
          }
        }
      }

      if(k-burn <=0 & ((burn-k)/n.lag)%%50 == 0){

        ii <- (burn-k)/n.lag; print(paste("burn-in iteration",ii))

        }

      #save values for each iteration
      if(k-burn > 0){
        if(((k-burn)%%n.lag)==0){

          ii <- (k-burn)/n.lag
          if((ii%%50)==0) {print(paste("iteration",ii))}

          Beta.I[ii,]  <- beta.I
          Beta.II[ii,] <- beta.II
          B.I[,,ii]    <- b.I
          B.II[,,ii]   <- b.II
          VCovI[,,ii]  <- solve(OmegaI)
          VCovII[,,ii] <- solve(OmegaII)
          Y1[[ii]] <- YI
          Y2[[ii]] <- YII
          circular.ri[,ii] <- int_rand_c

        }
      }
    }

  }else if(q1 == 1 & q2 >= 1){

    ####real iterations####
    for(k in 1:(tm))
    {

      for(i in 1:mm$N){
        YI[[i]]  <- mm$R[[i]]*cos(mm$theta[[i]])
        YII[[i]] <- mm$R[[i]]*sin(mm$theta[[i]])
      }

      #simulations for beta, bi's and Omega

      beta.I    <- betaBlock.fRI(OmegaI, YI, mm$XI, mm$ZI, p1, A1, mm$N)
      beta.II   <- betaBlock.fRS(OmegaII, YII, mm$XII, mm$ZII, p2, A2, mm$N)
      b.I       <- b.f(OmegaI, beta.I, YI, mm$XI, q1, mm$ZI, mm$ZtZI, mm$N)
      b.II      <- b.f(OmegaII, beta.II, YII, mm$XII, q2, mm$ZII, mm$ZtZII, mm$N)
      OmegaI    <- Omega.f(b.I, B1, v1, q1, mm$N)
      OmegaII   <- Omega.f(b.II ,B2, v2, q2, mm$N)

      for(i in 1:mm$N){

        if(!("(Intercept)" %in% colnames(mm$XI[[i]])) & !("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i,1], 0+b.I[i])

        }else if(!("(Intercept)" %in% colnames(mm$XI[[i]]))){

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i,1], 0+b.I[i])

        }else if(!("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i,1], beta.I[1]+b.I[i])

        }else{

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i,1], beta.I[1]+b.I[i])

        }

        for(j in 1:mm$no.Meas[i]){

          t.aux     <- mm$theta[[i]][j]
          mu.ij.I   <- c(beta.I%*%mm$XI[[i]][j,] + b.I[i])
          mu.ij.II  <- c(beta.II%*%mm$XII[[i]][j,] + b.II[i,]%*%mm$ZII[[i]][j,])
          bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
          mm$R[[i]][j] <- slice_r_me(t.aux, mu.ij.I, mu.ij.II, mm$R[[i]][j])

          if(k-burn > 0){
            if(((k-burn)%%n.lag)==0){
              ii <- (k-burn)/n.lag

              norm2  <- mu.ij.I^2 + mu.ij.II^2
              c         <- (1 + ((bb*pnorm(bb))/dnorm(bb)))

              predictiva[[i]][ii,j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2))*c
            }
          }
        }
      }

      if(k-burn <=0 & ((burn-k)/n.lag)%%50 == 0){

        ii <- (burn-k)/n.lag; print(paste("burn-in iteration",ii))

        }

      #save values for each iteration
      if(k-burn > 0){
        if(((k-burn)%%n.lag)==0){

          ii <- (k-burn)/n.lag
          if((ii%%50)==0) {print(paste("iteration",ii))}

          Beta.I[ii,]  <- beta.I
          Beta.II[ii,] <- beta.II
          B.I[,,ii]    <- b.I
          B.II[,,ii]   <- b.II
          VCovI[,,ii]  <- solve(OmegaI)
          VCovII[,,ii] <- solve(OmegaII)
          Y1[[ii]] <- YI
          Y2[[ii]] <- YII
          circular.ri[,ii] <- int_rand_c

        }
      }
    }



  }else{

    ####real iterations####
    for(k in 1:(tm)){

      for(i in 1:mm$N){
        YI[[i]]  <- mm$R[[i]]*cos(mm$theta[[i]])
        YII[[i]] <- mm$R[[i]]*sin(mm$theta[[i]])
      }

      #simulations for beta, bi's and Omega

      beta.I    <- betaBlock.fRS(OmegaI, YI, mm$XI, mm$ZI, p1, A1, mm$N)
      beta.II   <- betaBlock.fRS(OmegaII, YII, mm$XII, mm$ZII, p2, A2, mm$N)
      b.I       <- b.f(OmegaI, beta.I, YI, mm$XI, q1, mm$ZI, mm$ZtZI, mm$N)
      b.II      <- b.f(OmegaII, beta.II, YII, mm$XII, q2, mm$ZII, mm$ZtZII, mm$N)
      OmegaI    <- Omega.f(b.I, B1, v1, q1, mm$N)
      OmegaII   <- Omega.f(b.II ,B2, v2, q2, mm$N)

      for(i in 1:mm$N){

        if(!("(Intercept)" %in% colnames(mm$XI[[i]])) & !("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i,1], 0+b.I[i,1])

        }else if(!("(Intercept)" %in% colnames(mm$XI[[i]]))){

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i,1], 0+b.I[i,1])

        }else if(!("(Intercept)" %in% colnames(mm$XII[[i]]))){

          int_rand_c[i] <- atan2(0+b.II[i,1], beta.I[1]+b.I[i,1])

        }else{

          int_rand_c[i] <- atan2(beta.II[1]+b.II[i,1], beta.I[1]+b.I[i,1])

        }

        for(j in 1:mm$no.Meas[i]){

          t.aux     <- mm$theta[[i]][j]
          mu.ij.I   <- c(beta.I%*%mm$XI[[i]][j,] + b.I[i,]%*%mm$ZI[[i]][j,])
          mu.ij.II  <- c(beta.II%*%mm$XII[[i]][j,] + b.II[i,]%*%mm$ZII[[i]][j,])
          bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
          mm$R[[i]][j] <- slice_r_me(t.aux, mu.ij.I, mu.ij.II, mm$R[[i]][j])

          if(k-burn > 0){
            if(((k-burn)%%n.lag)==0){
              ii <- (k-burn)/n.lag

              norm2  <- mu.ij.I^2 + mu.ij.II^2
              c         <- (1 + ((bb*pnorm(bb))/dnorm(bb)))

              predictiva[[i]][ii,j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2))*c
            }
          }
        }
      }


      if(k-burn <=0 & ((burn-k)/n.lag)%%50 == 0){

        ii <- (burn-k)/n.lag; print(paste("burn-in iteration",ii))

        }

      #save values for each iteration
      if(k-burn > 0){
        if(((k-burn)%%n.lag)==0){

          ii <- (k-burn)/n.lag
          if((ii%%50)==0) {print(paste("iteration",ii))}

          Beta.I[ii,]  <- beta.I
          Beta.II[ii,] <- beta.II
          B.I[,,ii]    <- b.I
          B.II[,,ii]   <- b.II
          VCovI[,,ii]  <- solve(OmegaI)
          VCovII[,,ii] <- solve(OmegaII)
          Y1[[ii]] <- YI
          Y2[[ii]] <- YII
          circular.ri[,ii] <- int_rand_c

        }
      }
    }
  }

  output <- list(Beta.I = Beta.I, Beta.II = Beta.II, B.I = B.I, B.II = B.II,
                 VCovI = VCovI, VCovII = VCovII, predictiva = predictiva,
                 circular.ri = circular.ri, N = mm$N,
                 its = its, n.lag = n.lag, burn = burn,
                 p1 = p1, p2 = p2, q1 = q1, q2 = q2)

  summary.stats <- summe(output, mm)

  class(output) <- c("bpnme", class(output))

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

