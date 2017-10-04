#' Traceplots
#'
#' Traceplot function for a \code{bpnr object} or \code{bpnme object}.
#'
#' @param object an object used to select a method
#' @param parameter one of c("B1", "B2", "a.x", "a.c", "b.c", "SAM", "AS", "SSDO",
#'   "circ.diff") to indicate for which parameter a traceplot is required
#' @param variable a character string with variable name(s) to indicate for
#'   which variable(s) a traceplot is required
#'
#' @export
#'
#'

traceplot <- function(object, parameter = "SAM", variable){

UseMethod("traceplot", object)

}

#' Bayes Factors
#'
#' \code{BF} gives bayes factors for inequality constrained hypotheses on
#' circular mean differences
#'
#' @param object an object used to select a method
#' @param hypothesis the inequality constrained hypothesis to test
#' @param type type of hypothesis to test c("anchor", "isotropic")
#' @param ... further arguments passed to or from other methods
#'
#' @details the methods \link[pnreg]{BFc.bpnr} and
#'   \link[pnreg]{BFc.bpnme} have their own help page
#'
#' @export
#'

BFc <- function(object, hypothesis, type = "anchor"){

  UseMethod("BFc", object)

}

#' Model fit
#'
#' \code{fit} gives several model fit statistics
#'
#' @param object an object used to select a method
#'
#' @details the methods \link[pnreg]{fit.bpnr} and
#'   \link[pnreg]{fit.bpnme} have their own help page
#'
#' @export
#'

fit <- function(object){

  UseMethod("fit", object)

}

#' Linear coefficients
#'
#' \code{coef_lin} gives posterior summaries of the linear coefficients
#'
#' @param object an object used to select a method
#'
#' @details the methods \link[pnreg]{coef_lin.bpnr} and
#'   \link[pnreg]{coef_lin.bpnme} have their own help page
#'
#' @export
#'

coef_lin <- function(object){

  UseMethod("coef_lin", object)

}

#' Circular coefficients
#'
#' \code{coef_circ} gives posterior summaries of the circular coefficients
#'
#' @param object an object used to select a method
#' @param type one of c("continuous", "categorical") to get either the
#'   coefficients for the continuous or categorical predictor variables
#'
#' @details the methods \link[pnreg]{coef_circ.bpnr} and
#'   \link[pnreg]{coef_circ.bpnme} have their own help page
#'
#'
#' @export
#'

coef_circ <- function(object, type = "continuous"){

  UseMethod("coef_circ", object)

}

#' Random effect variances
#'
#' \code{coef_ran} gives posterior summaries of the circular or linear random effect variances
#'
#' @param object an object used to select a method
#' @param type one of c("linear", "circular") to get either the
#'   linear or circular variances of the random effects
#'
#' @details the method \link[pnreg]{coef_ran.bpnme} has its own help page
#'
#'
#' @export
#'

coef_ran <- function(object, type = "linear"){

  UseMethod("coef_ran", object)

}

#' Bayes Factors for a Bayesian circular regression model
#'
#' Outputs Bayes Factors for the circular differences between several levels of
#' a categorical variable and the baseline
#'
#' @param object a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#' @param hypothesis the inequality constrained hypothesis to test
#' @param type type of hypothesis to test c("anchor", "isotropic")
#'
#' @return Bayes Factors for inequality constrained hypothesis on mean differences
#'
#' @method BFc bpnr
#'
#' @export
#'

BFc.bpnr <- function(object, hypothesis, type = "anchor"){

  if(type == "anchor"){

    hypothesis2 <- gsub("[[:space:]]", "", hypothesis)
    variables <- unlist(strsplit(hypothesis2, split = "[<&>]"))

    all <- unlist(strsplit(hypothesis, split = "[[:space:]]"))
    logical.index <- attr(regexpr("[<&>]", all), "match.length")
    logical <- all[logical.index == 1]

    no.diff <- length(variables)

    data <- as.data.frame(object$circ.diff)[,variables]

    colnames(data) <- letters[1:no.diff]

    new.hyp <- c()

    for(i in 1:no.diff){

      name <- letters[i]

      if(i <= no.diff-1){

        new.hyp <- c(new.hyp, name, logical[i])

      }else{

        new.hyp <- c(new.hyp, name)

      }

    }

    new.hyp <- paste(new.hyp, collapse = "")


    subset <- subset(data, eval(parse(text = new.hyp)))

    fit <- nrow(subset)/nrow(data)

    complexity <- 1/factorial(length(unique(variables)))

    BF <- fit/complexity

    out <- matrix(NA, 2L, 3L)
    rownames(out) <- c("hypothesis", "complement")
    colnames(out) <- c("fit", "complexity", "Bayes Factor")

    out[1,] <- c(fit, complexity, BF)
    out[2,] <- c(1-fit, 1-complexity, (1-fit)/(1-complexity))


  }else if(type == "isotropic"){



  }

  return(out)

}

#' Bayes Factors for a Bayesian circular mixed-effects model
#'
#' Outputs Bayes Factors for inequality constrained hypotheses on the circular
#' differences between several levels of a categorical variable and the baseline
#'
#' @param object a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#' @param hypothesis the inequality constrained hypothesis to test
#' @param type type of hypothesis to test c("anchor", "isotropic")
#'
#' @return Bayes Factors for inequality constrained hypothesis on mean differences
#'
#' @method BFc bpnme
#'
#' @export
#'

BFc.bpnme <- function(object, hypothesis, type = "anchor"){

  if(type == "anchor"){

    hypothesis2 <- gsub("[[:space:]]", "", hypothesis)
    variables <- unlist(strsplit(hypothesis2, split = "[<&>]"))

    all <- unlist(strsplit(hypothesis, split = "[[:space:]]"))
    logical.index <- attr(regexpr("[<&>]", all), "match.length")
    logical <- all[logical.index == 1]

    no.diff <- length(variables)

    data <- as.data.frame(object$circ.diff)[,variables]

    colnames(data) <- letters[1:no.diff]

    new.hyp <- c()

    for(i in 1:no.diff){

      name <- letters[i]

      if(i <= no.diff-1){

        new.hyp <- c(new.hyp, name, logical[i])

      }else{

        new.hyp <- c(new.hyp, name)

      }

    }

    new.hyp <- paste(new.hyp, collapse = "")


    subset <- subset(data, eval(parse(text = new.hyp)))

    fit <- nrow(subset)/nrow(data)

    complexity <- 1/factorial(length(unique(variables)))

    BF <- fit/complexity

    out <- matrix(NA, 2L, 3L)
    rownames(out) <- c("hypothesis", "complement")
    colnames(out) <- c("fit", "complexity", "Bayes Factor")

    out[1,] <- c(fit, complexity, BF)
    out[2,] <- c(1-fit, 1-complexity, (1-fit)/(1-complexity))


  }else if(type == "isotropic"){



  }

  return(out)

}

#' Predicted values for a Bayesian circular regression model
#'
#' Outputs predicted values for a Bayesian circular regression model for each iteration
#' of the MCMC sampler
#'
#' @param object a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#' @param ... further arguments passed to or from other methods
#'
#' @return a matrix (rows = N, columns = iterations) containing predicted values
#'   for the circular outcome for each iteration of the MCMC sampler.
#'
#' @method predict bpnr
#'
#' @export
#'

predict.bpnr <- function(object, ...){

  YI <- object$mm$XI%*%t(object$B1)
  YII <- object$mm$XII%*%t(object$B2)

  theta <- atan2(YII, YI)%%(2*pi)

  return(theta)

}

#' Predicted values for a Bayesian circular mixed-effects model
#'
#' Outputs predicted values for a Bayesian circular mixed-effects model for each iteration
#' of the MCMC sampler
#'
#' @param object a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#' @param ... further arguments passed to or from other methods
#'
#' @return a matrix (rows = N, columns = iterations) containing predicted values
#'   for the circular outcome for each iteration of the MCMC sampler.
#'
#' @method predict bpnme
#'
#' @export
#'

predict.bpnme <- function(object, ...){

  YI <- list()
  YII <- list()
  theta <- list()

  for(i in 1:object$mm$N){

    YI[[i]] <- object$mm$XI[[i]] %*% t(object$Beta.I) + object$mm$ZI[[i]]%*%object$B.I[i,,]
    YII[[i]] <- object$mm$XII[[i]]%*%t(object$Beta.II) + object$mm$ZII[[i]]%*%object$B.II[i,,]

    theta[[i]] <- atan2(YII[[i]], YI[[i]])

  }

  YI <- do.call(rbind, YI)
  YII <- do.call(rbind, YII)
  theta <- do.call(rbind, theta)%%(2*pi)

  return(theta)

}

#' Residuals for a Bayesian circular regression model
#'
#' Outputs residuals for a Bayesian circular regression model for each iteration
#' of the MCMC sampler
#'
#' @param object a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#' @param type the type of residuals, one of \code{c("arc", "cos")}. The "arc"
#'   residuals are based on a computation of the circular arc length between
#'   predicted value and original outcome. The "cos" residuals are based on the
#'   cosine of the difference between predicted value and original outcome.
#' @param ... further arguments passed to or from other methods
#'
#' @return a matrix (rows = N, columns = iterations) containing residuals for
#'   each iteration of the MCMC sampler.
#'
#' @method residuals bpnr
#'
#' @export
#'

residuals.bpnr <- function(object, type = "arc", ...){

  pred.val <- predict(object)
  diff <- pred.val - as.numeric(object$theta)
  sign <- sign(sin(diff))

  if(type == "arc"){

    residuals <- (pi - abs(pi - abs(diff)))*sign

  }else if(type == "cos"){

    residuals <- cos(abs(diff))

  }

  return(residuals)

}

#' Residuals for a Bayesian circular mixed-effects model
#'
#' Outputs residuals for a Bayesian circular mixed-effects model for each iteration
#' of the MCMC sampler
#'
#' @param object a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#' @param type the type of residuals, one of \code{c("arc", "cos")}. The "arc"
#'   residuals are based on a computation of the circular arc length between
#'   predicted value and original outcome. The "cos" residuals are based on the
#'   cosine of the difference between predicted value and original outcome.
#' @param ... further arguments passed to or from other methods
#'
#' @return a matrix (rows = N, columns = iterations) containing residuals for
#'   each iteration of the MCMC sampler.
#'
#' @method residuals bpnme
#'
#' @export
#'

residuals.bpnme <- function(object, type = "arc", ...){

  pred.val <- predict(object)

  theta <- do.call(rbind, object$mm$theta)

  diff <- pred.val - as.numeric(theta)
  sign <- sign(sin(diff))

  if(type == "arc"){

    residuals <- (pi - abs(pi - abs(diff)))*sign

  }else if(type == "cos"){

    residuals <- cos(abs(diff))

  }

  return(residuals)

}

#' Obtain random effect variances of a Bayesian circular mixed-effects model
#'
#' Gives posterior summaries of the circular or linear random effect variances
#'
#' @param object a \code{bpnme object} obtained from the function \code{bpnme()}
#' @param type one of c("linear", "circular") to get either the
#'   linear or circular variances of the random effects
#'
#' @return A matrix with posterior summaries of the random effect variances
#'
#'
#' @export
#'

coef_ran.bpnme <- function(object, type = "linear"){

  if(type == "linear"){

    return(rbind(object$lin.res.varrand.I, object$lin.res.varrand.II))

  }else if(type == "circular"){

    return(object$circ.res.varrand)

  }else{

    stop("type not recognized")

  }

}

#' Obtain the linear coefficients of a Bayesian circular regression model
#'
#' Gives the coefficients tables of the linear coefficients for a circular
#' regression model
#'
#' @param object a \code{bpnr object} obtained from the function \code{bpnr()}
#'
#' @return A matrix with posterior summaries of the linear coefficients in a
#'   Bayesian circular regression model
#'
#' @method coef_lin bpnr
#'
#' @export
#'
#'

coef_lin.bpnr <- function(object){

  return(rbind(object$lin.coef.I, object$lin.coef.II))

}

#' Obtain the linear coefficients of a Bayesian circular mixed-effects model
#'
#' Gives the coefficients tables of the linear coefficients for a Bayesian
#' circular mixed-effects model
#'
#' @param object a \code{bpnme object} obtained from the function \code{bpnme()}
#'
#' @return A matrix with posterior summaries of the linear coefficients in a
#'   Bayesian circular mixed-effects model
#'
#' @method coef_lin bpnme
#'
#' @export
#'
#'

coef_lin.bpnme <- function(object){

  return(rbind(object$lin.coef.I, object$lin.coef.II))

}

#' Obtain the circular coefficients of a Bayesian circular regression model
#'
#' Gives the coefficients tables of the circular coefficients for a Bayesian circular regression model
#'
#' @param object a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#' @param type one of c("continuous", "categorical") to get either the
#'   coefficients for the continuous or categorical predictor variables
#'
#' @return A matrix with posterior summaries of the circular coefficients in a
#'   Bayesian circular regression model
#'
#' @method coef_circ bpnr
#'
#' @export
#'
#'

coef_circ.bpnr <- function(object, type = "continuous"){

  if(type == "continuous"){

    if(class(object$circ.coef) == "character"){

      return(object$circ.coef)

    }else{



      a.x = object$circ.coef[,1:5]
      a.c = object$circ.coef[,6:10]
      b.c = object$circ.coef[,11:15]
      AS = object$circ.coef[,16:20]
      SAM = object$circ.coef[,21:25]
      SSDO = object$circ.coef[,26:30]

      coefficients <- array(c(a.x, a.c, b.c, AS, SAM, SSDO),
                            dim = c(nrow(object$circ.coef), 5, 6))
      dimnames(coefficients) <- list(rownames(object$circ.coef),
                                     c("mean", "mode", "sd", "LB HPD", "UB HPD"),
                                     c("ax", "ac", "bc", "AS", "SAM", "SSDO"))

      return(coefficients)

    }


  }else if(type == "categorical"){

    return(object$circ.coef.cat)

  }else{

    stop("type not recognized")

  }

}

#' Obtain the circular coefficients of a Bayesian circular mixed-effects model
#'
#' Gives the coefficients tables of the circular coefficients for a Bayesian circular mixed-effects model
#'
#' @param object a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#' @param type one of c("continuous", "categorical") to get either the
#'   coefficients for the continuous or categorical predictor variables
#'
#' @return A matrix with posterior summaries of the circular coefficients in a
#'   Bayesian circular mixed-effects model
#'
#' @method coef_circ bpnme
#'
#' @export
#'
#'

coef_circ.bpnme <- function(object, type = "continuous"){

  if(type == "continuous"){

    if(class(object$circ.coef) == "character"){

      return(object$circ.coef)

    }else{



      a.x = object$circ.coef[,1:5]
      a.c = object$circ.coef[,6:10]
      b.c = object$circ.coef[,11:15]
      AS = object$circ.coef[,16:20]
      SAM = object$circ.coef[,21:25]
      SSDO = object$circ.coef[,26:30]

      coefficients <- array(c(a.x, a.c, b.c, AS, SAM, SSDO),
                            dim = c(nrow(object$circ.coef), 5, 6))
      dimnames(coefficients) <- list(rownames(object$circ.coef),
                                     c("mean", "mode", "sd", "LB HPD", "UB HPD"),
                                     c("ax", "ac", "bc", "AS", "SAM", "SSDO"))

      return(coefficients)

    }


  }else if(type == "categorical"){

    return(object$circ.coef.cat)

  }else{

    stop("type not recognized")

  }

}

#' Model fit for a Bayesian circular regression model
#'
#' Outputs several model fit statistics for the Bayesian circular regression model
#'
#' @param object a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#'
#' @return a matrix containing the computed log pointwise predictive density
#'   (lppd), Deviance Information Criterion (DIC), an alternative version of the
#'   DIC (DIC_alt), and the Watanabe-Akaike Information Criterion computed in
#'   two different ways (WAIC, WAIC2). The matrix also contains the number of
#'   parameters or 'effective number' of parameters that the several statistics
#'   are based on. Computation of the criteria is done accoring to Gelman et.al
#'   (2014) in *Bayesian Data Analysis*.
#'
#' @method fit bpnr
#'
#' @export
#'

fit.bpnr <- function(object){

  mf <- matrix(NA, 5L, 2L)
  colnames(mf) <- c("Statistic", "Parameters")
  rownames(mf) <- c("lppd", "DIC", "DIC.alt", "WAIC", "WAIC2")
  mf[,1] <- unlist(object$model.fit[c("lppd", "DIC", "DIC_alt", "WAIC", "WAIC2")])
  mf[,2] <- c(object$p1 + object$p2,
              unlist(object$model.fit[c("pD", "pV", "pWAIC", "pWAIC2")]))
  as.data.frame(mf)

}

#' Model fit for a Bayesian circular mixed-effects model
#'
#' Outputs several model fit statistics for the Bayesian circular mixed-effects model
#'
#' @param object a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#'
#' @return a matrix containing the computed log pointwise predictive density
#'   (lppd), Deviance Information Criterion (DIC), an alternative version of the
#'   DIC (DIC_alt), and the Watanabe-Akaike Information Criterion computed in
#'   two different ways (WAIC, WAIC2). The matrix also contains the number of
#'   parameters or 'effective number' of parameters that the several statistics
#'   are based on. Computation of the criteria is done accoring to Gelman et.al
#'   (2014) in *Bayesian Data Analysis*.
#'
#' @method fit bpnme
#'
#' @export
#'

fit.bpnme <- function(object){

  mf <- matrix(NA, 5L, 2L)
  colnames(mf) <- c("Statistic", "Parameters")
  rownames(mf) <- c("lppd", "DIC", "DIC.alt", "WAIC", "WAIC2")
  mf[,1] <- object$model.fit[dimnames(object$model.fit)[[2]] %in% c("lppd", "DIC", "DICalt", "WAIC", "WAIC2")]
  mf[,2] <- c(object$p1 + object$p2,
              object$model.fit[dimnames(object$model.fit)[[2]] %in% c("pD", "pV", "pWAIC", "pWAIC2")])
  as.data.frame(mf)

}


#' Print output from a Bayesian circular regression model
#'
#' Prints selected output from a Bayesian circular regression model
#'
#' @param x a \code{bpnr object} obtained from the function
#'   \code{bpnr()}
#' @param ... further arguments passed to or from other methods
#'
#' @return a print of selected output from a Bayesian circular regression model
#'
#' @method print bpnr
#'
#' @export
#'

print.bpnr <- function(x, ...){

  cat("Projected Normal Regression \n\n")
  cat("Model \n\n")

  cat("Call: \n", paste(deparse(x$Call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("MCMC: \n", paste("iterations = ", x$its, "\n",
                        "burn-in = ", x$burn, "\n",
                        "lag = ", x$lag,
                        sep = ""),
      "\n\n", sep = "")

  cat("Model Fit: \n")
  mf <- matrix(NA, 5L, 2L)
  colnames(mf) <- c("Statistic", "Parameters")
  rownames(mf) <- c("lppd", "DIC", "DIC.alt", "WAIC", "WAIC2")
  mf[,1] <- unlist(x$model.fit[c("lppd", "DIC", "DIC_alt", "WAIC", "WAIC2")])
  mf[,2] <- c(x$p1 + x$p2,
              unlist(x$model.fit[c("pD", "pV", "pWAIC", "pWAIC2")]))
  print(mf)
  cat("\n\n")

  cat("Linear Coefficients \n\n")

  cat("Component I: \n")
  print(x$lin.coef.I)
  cat("\n")

  cat("Component II: \n")
  print(x$lin.coef.II)
  cat("\n\n")

  cat("Circular Coefficients \n\n")

  cat("Continuous variables: \n")
  if(class(x$circ.coef) == "character"){
    print(x$circ.coef)
  }else{
    print(x$circ.coef[,1:5])
    cat("\n")
    print(x$circ.coef[,6:10])
    cat("\n")
    print(x$circ.coef[,11:15])
    cat("\n")
    print(x$circ.coef[,16:20])
    cat("\n")
    print(x$circ.coef[,21:25])
    cat("\n")
    print(x$circ.coef[,26:30])
  }
  cat("\n")

  cat("Categorical variables: \n")
  print(x$circ.coef.cat)
  cat("\n")

}

#' Print output from a Bayesian circular mixed-effects model
#'
#' Prints selected output from a Bayesian circular mixed-effects model
#'
#' @param x a \code{bpnme object} obtained from the function
#'   \code{bpnme()}
#' @param ... further arguments passed to or from other methods
#'
#' @return a print of selected output from a Bayesian circular mixed-effects model
#'
#' @method print bpnme
#'
#' @export
#'

print.bpnme <- function(x, ...){

  cat("Projected Normal Mixed Effects \n\n")
  cat("Model \n\n")

  cat("Call: \n", paste(deparse(x$Call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("MCMC: \n", paste("iterations = ", x$its, "\n",
                        "burn-in = ", x$burn, "\n",
                        "lag = ", x$n.lag,
                        sep = ""),
      "\n\n", sep = "")

  cat("Model Fit: \n")
  mf <- matrix(NA, 5L, 2L)
  colnames(mf) <- c("Statistic", "Parameters")
  rownames(mf) <- c("lppd", "DIC", "DIC.alt", "WAIC", "WAIC2")
  mf[,1] <- x$model.fit[dimnames(x$model.fit)[[2]] %in% c("lppd", "DIC", "DICalt", "WAIC", "WAIC2")]
  mf[,2] <- c(x$p1 + x$p2,
              x$model.fit[dimnames(x$model.fit)[[2]] %in% c("pD", "pV", "pWAIC", "pWAIC2")])
  print(mf)
  cat("\n\n")

  cat("Fixed Effects \n\n")

  cat("Linear Coefficients \n\n")

  cat("Component I: \n")
  print(x$lin.coef.I)
  cat("\n")

  cat("Component II: \n")
  print(x$lin.coef.II)
  cat("\n\n")

  cat("Circular Coefficients \n\n")

  cat("Continuous variables: \n")
  if(class(x$circ.coef) == "character"){
    print(x$circ.coef)
  }else{
    print(x$circ.coef[,1:5])
    cat("\n")
    print(x$circ.coef[,6:10])
    cat("\n")
    print(x$circ.coef[,11:15])
    cat("\n")
    print(x$circ.coef[,16:20])
    cat("\n")
    print(x$circ.coef[,21:25])
    cat("\n")
    print(x$circ.coef[,26:30])
  }
  cat("\n")

  cat("Categorical variables: \n")
  print(x$circ.coef.cat)
  cat("\n\n")

  cat("Random Effects \n\n")

  cat("Linear Coefficients \n\n")

  cat("Component I: \n")
  print(x$lin.res.varrand.I)
  cat("\n")

  cat("Component II: \n")
  print(x$lin.res.varrand.II)
  cat("\n\n")

  cat("Circular Coefficients \n\n")
  print(x$circ.res.varrand)

}

#' Traceplots for a Bayesian circular regression model
#'
#' General plot function for a \code{bpnr object}.
#'
#' @param object a \code{bpnr object} obtained from the function \code{bpnr()}
#' @param parameter one of c("B1", "B2", "a.x", "a.c", "b.c", "SAM", "AS", "SSDO",
#'   "circ.diff") to indicate for which parameter a traceplot is required
#' @param variable a character string with variable name(s) to indicate for
#'   which variable(s) a traceplot is required
#'
#' @method traceplot bpnr
#'
#' @export
#'
#'

traceplot.bpnr <- function(object, type = "trace", parameter = "SAM", variable = NULL){

 if(is.null(variable)){

   text <- paste0(as.character(bquote(object)), "$", parameter)
   plot.ts(eval(parse(text = text)), xlab = "Iteration", main = "Traceplots")

 }else{

   text <- paste0(as.character(bquote(object)), "$", parameter)
   plot.ts(eval(parse(text = text))[,variable], ylab = variable, xlab = "Iteration", main = "Traceplots")

 }



}

#' Traceplots for a Bayesian circular mixed-effects model
#'
#' General plot function for a \code{bpnme object}.
#'
#' @param object a \code{bpnme object} obtained from the function \code{bpnme()}
#' @param parameter one of c("Beta.I", "Beta.II", a.x", "a.c", "b.c", "SAM",
#'   "AS", "SSDO", "circ.diff") to indicate for which parameter a traceplot is
#'   required
#' @param variable if effect = "fixed", a character string with variable name(s)
#'   to indicate for which variable(s) a traceplot is required
#'
#' @method traceplot bpnme
#'
#' @export
#'
#'

traceplot.bpnme <- function(object, parameter = "SAM", variable){

  if(effect == "fixed"){

    text <- paste0(as.character(bquote(object)), "$", parameter)
    plot.ts(eval(parse(text = text))[,variable], ylab = variable, xlab = "Iteration", main = "Traceplots")

  }

  if(effect == "random"){


  }

}
