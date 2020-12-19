#' Compute utmu
#'
#' @param t current outcome value
#' @param mu1 current predicted linear mean of the first component
#' @param mu2 current predicted linear mean of the second component
#'
#' @keywords internal
#'

Dbd <- function(t, mu1, mu2){ cos(t)*mu1 + sin(t)*mu2 }

#' Check whether a variable is categorical
#'
#' \code{cat_check} checks whether a vector contains only 0 and 1 and thus is a
#' dummy variable.
#'
#' @param x vector containing values of a variable.
#'
#' @return TRUE, if the vector only contains 0 and 1, FALSE is the vector
#'   contains other values.
#'
#' @keywords internal
#'

cat_check <- function(x){
  length(which(!x%in%c(0,1))) == 0
}

#' Compute circular coefficients from linear coefficients
#'
#' \code{circ_coef} computes the coordinates of the inflection point of
#' a circular effect, the slope at the inflection point and the unsigned and
#' signed shortest distance to the origin.
#'
#' @param a1 intercept of the first linear component.
#' @param a2 intercept of the second linear component.
#' @param b1 slope of the first linear component.
#' @param b2 slope of the second linear component.
#'
#' @return A dataframe containing the coordinates of the inflection point of a
#'   circular effect, the slope at the inflection point and the unsigned and
#'   signed shortest distance to the origin.
#'

circ_coef <- function(a1, a2, b1, b2){

  ax <- - (a1*b1 + a2*b2) / (b1^2 + b2^2)
  ac <- atan2(a2 + b2*ax, a1 + b1*ax)
  bc <- (tan(atan2(a2, a1)) - ac) / -ax

  SDO <- sqrt( (a1 + b1)^2 + (a2 + b2)^2)
  SSDO <- sign(sin(ac - atan2(b2, b1)))*SDO

  Output <- cbind(ax, ac, bc, SDO, SSDO)
  colnames(Output) <- c("ax", "ac", "bc", "SDO", "SSDO")

  return(as.data.frame(Output))
}

#' Compute the standard deviation of a vector of circular data
#'
#' @param theta a circular variable in radians or degrees.
#' @param units measurement units of the circular variable c("radians",
#'   "degrees").
#'
#' @examples
#' library(bpnreg)
#' sd_circ(subset(Motor, Cond == "exp")$PhaseDiff, units = "degrees")
#'
#' @export

sd_circ <- function(theta, units = "radians"){

  if(!(units == "degrees" | units == "radians")){

    stop("Invalid units argument")

  }

  if(units == "radians"){

    sqrt(-2*log(rho_circ(theta)))

  }else if(units == "degrees"){

    sqrt(-2*log(rho_circ(theta, units = "degrees")))*(180/pi)

  }


}

#' Compute the mean of a vector of circular data
#'
#' @inheritParams sd_circ
#'
#' @examples
#' library(bpnreg)
#' mean_circ(subset(Motor, Cond == "exp")$PhaseDiff, units = "degrees")
#'
#' @export

mean_circ <- function(theta, units = "radians"){

  if(!(units == "degrees" | units == "radians")){

    stop("Invalid units argument")

  }

  if(units == "radians"){

    theta_bar(theta)

  }else if(units == "degrees"){

    theta_bar((theta*(pi/180)))*(180/pi)

  }

}

#' Compute the mean resultant length of a vector of circular data
#'
#' @inheritParams sd_circ
#'
#' @examples
#' library(bpnreg)
#' rho_circ(subset(Motor, Cond == "exp")$PhaseDiff, units = "degrees")
#'
#' @export

rho_circ <- function(theta, units = "radians"){

  if(!(units == "degrees" | units == "radians")){

    stop("Invalid units argument")

  }

  if(units == "radians"){

    rho(theta)$rho

  }else if(units == "degrees"){

    rho(theta*(pi/180))$rho

  }

}

#' Compute the mode of a vector of linear data
#'
#' @param x a vector of linear data
#'
#' @examples
#' library(bpnreg)
#' mode_est(Motor$AvAmp)
#'
#' @export

mode_est <- function(x){hmode(x, 0.1)}

#' Compute the 95 percent HPD of a vector of linear data
#'
#' @inheritParams mode_est
#'
#' @examples
#' library(bpnreg)
#' hpd_est(Motor$AvAmp)
#'
#' @export

hpd_est <- function(x){hmodeci(x, 0.95)}

#' Compute the mode of a vector of circular data
#'
#' @param x a vector of circular data
#'
#' @examples
#' library(bpnreg)
#' mode_est_circ(subset(Motor, Cond = "exp")$Phaserad)
#'
#' @export

mode_est_circ <- function(x){hmodeC(x, 0.1)}

#' Compute the 95 percent HPD of a vector of circular data
#'
#' @examples
#' library(bpnreg)
#' hpd_est_circ(subset(Motor, Cond = "exp")$Phaserad)
#'
#' @inheritParams mode_est_circ
#'
#' @export

hpd_est_circ <- function(x){hmodeciC(x, 0.95)}

#' Create model matrices circular regression
#'
#' @param pred.I model equation for effects of component 1
#' @param data the dataframe used for analysis
#' @param pred.II model equation for effects of component 2
#'
#'

mmr <- function(pred.I, data, pred.II){

  theta <- model.frame(pred.I, data)[,1]

  #Warning message if theta is not measures in radians:

  if(!(max(unlist(theta)) < 2*pi & min(unlist(theta)) >= 0) & !(max(unlist(theta)) < pi & min(unlist(theta)) >= -pi)){
    cat("The circular outcome may contain values that are out of range. \n")
    cat("Please check that the circular outcome is measured in radians.")
  }

  X_I   <- model.matrix(pred.I, data)
  X_II  <- model.matrix(pred.II, data)

  N <- nrow(X_I)

  R <- rep(1, N)

  list(XI = X_I, XII = X_II, N = N,
       theta = theta, R = R,
       pred.I = pred.I, pred.II = pred.II, data = data)

}

#' Create model matrices for a circular mixed-effects regression model
#'
#' @param pred.I model equation for effects of component 1
#' @param data the dataframe used for analysis
#' @param pred.II model equation for effects of component 2
#'
#'
#'

mmme <- function(pred.I, data, pred.II){

  fix_form.I <- pred.I
  fix_form.II <- pred.II
  ran_form.I <- pred.I
  ran_form.II <- pred.II

  nesting.I <- strsplit(sub("\\| ", "",
                            sub("^[^\\|]*", "",
                                attr(terms(reOnly(ran_form.I)), "term.labels"))), " ")[[1]]
  nesting.II <- strsplit(sub("\\| ", "",
                             sub("^[^\\|]*", "",
                                 attr(terms(reOnly(ran_form.II)), "term.labels"))), " ")[[1]]

  all_ran.I <- strsplit(attr(terms(reOnly(ran_form.I)), "term.labels"), " ")[[1]]
  all_ran.II <- strsplit(attr(terms(reOnly(ran_form.II)), "term.labels"), " ")[[1]]

  if(any(grepl("\\:", nesting.I)) | any(grepl("\\/", nesting.I)) | any(grepl("\\+", nesting.I)) |
     any(grepl("\\:", nesting.II)) | any(grepl("\\/", nesting.I)) | any(grepl("\\+", nesting.I))){
      stop("More than one nesting variable defined")
  }

  if(!all(sapply(data[, nesting.I], inherits, TRUE, what =  "numeric")) |
     !all(sapply(data[, nesting.II], inherits, TRUE, what =  "numeric"))){
    stop("Not all nesting variables are class numeric.")
  }

  for(i in nesting.I){
    if(length(which(all_ran.I == i)) > 1){
      stop("Nesting variables cannot be random slopes in the same model")
    }
  }

  for(i in nesting.II){
    if(length(which(all_ran.II == i)) > 1){
      stop("Nesting variables cannot be random slopes in the same model")
    }
  }

  RHSForm(fix_form.I) <- nobars(RHSForm(fix_form.I))
  RHSForm(fix_form.II) <- nobars(RHSForm(fix_form.II))
  RHSForm(ran_form.I) <- subbars(RHSForm(reOnly(ran_form.I)))
  RHSForm(ran_form.II) <- subbars(RHSForm(reOnly(ran_form.II)))

  lab_ran.I <- attr(terms(ran_form.I), "term.labels")
  lab_ran.II <- attr(terms(ran_form.II), "term.labels")
  lab_fix.I <- attr(terms(fix_form.I), "term.labels")
  lab_fix.II <- attr(terms(fix_form.II), "term.labels")

  if(length(which(lab_fix.I == nesting.I)) > 0 | length(which(lab_fix.II == nesting.II)) > 0){
    stop("Nesting variables cannot be fixed effects in the same model")
  }

  no_terms_ran.I <- length(lab_ran.I)
  no_terms_ran.II <- length(lab_ran.II)

  theta <- split(model.frame(ran_form.I, data)[,1], data[,lab_ran.I[no_terms_ran.I]])

  mm.I <- model.matrix(fix_form.I, data)
  mm.II <- model.matrix(fix_form.II, data)
  mm_ran.I <- model.matrix(ran_form.I, data)
  mm_ran.II <- model.matrix(ran_form.II, data)

  n_ran.I <- colnames(mm_ran.I)
  n_ran.II <- colnames(mm_ran.II)

  mm_ran.I <- as.matrix(mm_ran.I[,1:ncol(mm_ran.I)-1])
  colnames(mm_ran.I) <- n_ran.I[1:ncol(mm_ran.I)]

  mm_ran.II <- as.matrix(mm_ran.II[,1:ncol(mm_ran.II)-1])
  colnames(mm_ran.II) <- n_ran.II[1:ncol(mm_ran.II)]

  if(!"(Intercept)" %in% colnames(mm_ran.I) | !"(Intercept)" %in% colnames(mm_ran.II)){

    stop("No random intercept in the model")

  }

  if(!all(colnames(mm_ran.I) %in% colnames(mm.I)) |
     !all(colnames(mm_ran.II) %in% colnames(mm.II))){

    stop("Not all random effects have a corresponding fixed effect")

  }

  X_I   <- split(model.matrix(fix_form.I, data), data[,lab_ran.I[no_terms_ran.I]])
  X_II  <- split(model.matrix(fix_form.II, data), data[,lab_ran.II[no_terms_ran.II]])

  Z_I   <- split(model.matrix(ran_form.I, data)[,-(no_terms_ran.I + 1)], data[,lab_ran.I[no_terms_ran.I]])
  Z_II  <- split(model.matrix(ran_form.II, data)[, -(no_terms_ran.II + 1)], data[,lab_ran.II[no_terms_ran.II]])

  N <- nrow(unique(data[,lab_ran.I[no_terms_ran.I], drop = FALSE]))

  ZtZ_I   <- list()
  ZtZ_II  <- list()
  no.Meas <- c()
  R <- list()

  for(i in 1:N){

    no.Meas[i]    <- length(theta[[i]])
    R[[i]]        <- as.matrix(rep(1, length(theta[[i]])))
    theta[[i]]    <- as.matrix(theta[[i]])

    X_I[[i]]      <- matrix(X_I[[i]], no.Meas[i], length(X_I[[i]])/no.Meas[i], dimnames = list(NULL, colnames(mm.I)))
    X_II[[i]]     <- matrix(X_II[[i]], no.Meas[i], length(X_II[[i]])/no.Meas[i], dimnames = list(NULL, colnames(mm.II)))
    Z_I[[i]]      <- matrix(Z_I[[i]], no.Meas[i], length(Z_I[[i]])/no.Meas[i], dimnames = list(NULL, colnames(mm_ran.I)))
    Z_II[[i]]     <- matrix(Z_II[[i]], no.Meas[i], length(Z_II[[i]])/no.Meas[i], dimnames = list(NULL, colnames(mm_ran.II)))
    ZtZ_I[[i]]    <- t((Z_I[[i]]))%*%(Z_I[[i]])
    ZtZ_II[[i]]   <- t((Z_II[[i]]))%*%(Z_II[[i]])
  }

  #Warning message if theta is not measures in radians:

  if(!(max(unlist(theta)) < 2*pi & min(unlist(theta)) >= 0) & !(max(unlist(theta)) < pi & min(unlist(theta)) >= -pi)){
    cat("The circular outcome may contain values that are out of range. \n")
    cat("Please check that the circular outcome is measured in radians.")
  }

  list(XI = X_I, XII = X_II, ZI = Z_I, ZII = Z_II, N = N,
       no.Meas = no.Meas, ZtZI = ZtZ_I, ZtZII = ZtZ_II, theta = theta, R = R,
       pred.I = pred.I, pred.II = pred.II,
       fix_form.I = fix_form.I, fix_form.II = fix_form.II,
       ran_form.I = ran_form.I, ran_form.II = ran_form.II,
       mm.I = mm.I, mm.II = mm.II, mm_ran.I = mm_ran.I, mm_ran.II = mm_ran.II,
       data = data)

}

#' Compute summary and model fit statistics for the circular regression model
#'
#' @param output from the regression estimation function pnr()
#' @param mm output from the function mmr()
#'
#' @keywords internal
#'

sumr <- function(output, mm){

  var.names.I <-  colnames(mm$XI)
  var.names.II <-  colnames(mm$XII)

  if("(Intercept)" %in% colnames(mm$XI)){

    cat <- apply(mm$XI, 2, cat_check)
    var.cat.I <- colnames(mm$XI)[cat][-1]
    var.num.I <- colnames(mm$XI)[!cat]
    grandmeans.I <- colMeans(mm$XI)[-1]

  }else{

    cat <- apply(mm$XI, 2, cat_check)
    var.cat.I <- colnames(mm$XI)[cat]
    var.num.I <- colnames(mm$XI)[!cat]
    grandmeans.I <- colMeans(mm$XI)

  }

  if("(Intercept)" %in% colnames(mm$XII)){

    cat <- apply(mm$XII, 2, cat_check)
    var.cat.II <- colnames(mm$XII)[cat][-1]
    var.num.II <- colnames(mm$XII)[!cat]
    grandmeans.II <- colMeans(mm$XII)

  }else{

    cat <- apply(mm$XII, 2, cat_check)
    var.cat.II <- colnames(mm$XII)[cat]
    var.num.II <- colnames(mm$XII)[!cat]
    grandmeans.II <- colMeans(mm$XII)

  }

  var.num <- unique(c(var.num.I, var.num.II))
  var.cat <- unique(c(var.cat.I, var.cat.II))


  #Linear effects

  lin.res.I <- matrix(NA, ncol(output$beta1), 5)
  rownames(lin.res.I) <- c(var.names.I)
  colnames(lin.res.I) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")

  colnames(output$beta1) <- c(var.names.I)

  lin.res.II <- matrix(NA, ncol(output$beta2), 5)
  rownames(lin.res.II) <- c(var.names.II)
  colnames(lin.res.II) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")

  colnames(output$beta2) <- c(var.names.II)

  for(i in colnames(output$beta1)){

    lin.res.I[i,1] <- mean(output$beta1[,i])
    lin.res.I[i,2] <- mode_est(output$beta1[,i])
    lin.res.I[i,3] <- sd(output$beta1[,i])
    lin.res.I[i,4:5] <- hpd_est(output$beta1[,i])

  }

  for(i in colnames(output$beta2)){

    lin.res.II[i,1] <- mean(output$beta2[,i])
    lin.res.II[i,2] <- mode_est(output$beta2[,i])
    lin.res.II[i,3] <- sd(output$beta2[,i])
    lin.res.II[i,4:5] <- hpd_est(output$beta2[,i])

  }

  #Circular effects

  if(length(var.cat) >= 2){

    var.comb.cat <- combn(var.cat, 2)

    circ.res.means <- matrix(NA, length(var.cat) + ncol(var.comb.cat) + 1, 5)
    names <- sapply(seq_along(var.comb.cat[1,]),
                    function(w){paste(var.comb.cat[,w], sep = "", collapse = "")})
    rownames(circ.res.means) <- c("(Intercept)", var.cat, names)
    colnames(circ.res.means) <- c("mean", "mode", "sd", "LB", "UB")

    circ.diff <- matrix(NA, output$its, length(var.cat) + ncol(var.comb.cat))
    colnames(circ.diff) <- c(var.cat, names)

    circ.res.cat <- matrix(NA, length(var.cat) + ncol(var.comb.cat), 5)
    rownames(circ.res.cat) <- c(var.cat, names)
    colnames(circ.res.cat) <- c("mean", "mode", "sd", "LB", "UB")

  }else{

    var.comb.cat <- 0

    circ.res.means <- matrix(NA, length(var.cat) + 1, 5)
    rownames(circ.res.means) <- c("(Intercept)", var.cat)
    colnames(circ.res.means) <- c("mean", "mode", "sd", "LB", "UB")

    circ.diff <- matrix(NA, output$its, length(var.cat))
    colnames(circ.diff) <- var.cat

    circ.res.cat <- matrix(NA, length(var.cat), 5)
    rownames(circ.res.cat) <- var.cat
    colnames(circ.res.cat) <- c("mean", "mode", "sd", "LB", "UB")


  }


  circ.res <- matrix(NA, length(var.num), 5*6)
  rownames(circ.res) <- var.num
  colnames(circ.res) <- c("mean ax", "mode ax", "sd ax", "LB ax", "UB ax",
                          "mean ac", "mode ac", "sd ac", "LB ac", "UB ac",
                          "mean bc", "mode bc", "sd bc", "LB bc", "UB bc",
                          "mean AS", "mode AS", "sd AS", "LB AS", "UB AS",
                          "mean SAM", "mode SAM", "sd SAM", "LB SAM", "UB SAM",
                          "mean SSDO", "mode SSDO", "sd SSDO", "LB SSSO", "UB SSDO")

  if(length(var.cat) == 0){

    circ.res.cat <- "There are no categorical predictors in the model"
    circ.res.means <- "There are no categorical predictors in the model"
    circ.diff <- "There are no categorical predictors in the model"

  }else if(!("(Intercept)" %in% colnames(mm$XI)) & !("(Intercept)" %in% colnames(mm$XII))){

    circ.res.cat <- "There is no intercept in the model"
    circ.diff <- "There is no intercept in the model"

  }else if(all.equal(var.cat.I, var.cat.II) == TRUE){

    for(c in var.cat){

      if(!("(Intercept)" %in% colnames(mm$XI))){

        baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
        dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + output$beta1[, c])

      }else if(!("(Intercept)" %in% colnames(mm$XII))){

        baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
        dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

      }else{

        baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
        dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

      }

      diff <- baseline - dummy
      sign <- sign(sin(diff))
      circDiff <- (pi - abs(pi - abs(diff)))*sign

      circ.diff[,c] <- circDiff

      circ.res.cat[c,1] <- theta_bar(circDiff)
      circ.res.cat[c,2] <- mode_est_circ(circDiff)
      circ.res.cat[c,3] <- sd_circ(circDiff)
      circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

      circ.res.means[c,1] <- theta_bar(dummy)
      circ.res.means[c,2] <- mode_est_circ(dummy)
      circ.res.means[c,3] <- sd_circ(dummy)
      circ.res.means[c,4:5] <- hpd_est_circ(dummy)


    }

    circ.res.means[1,1] <- theta_bar(baseline)
    circ.res.means[1,2] <- mode_est_circ(baseline)
    circ.res.means[1,3] <- sd_circ(baseline)
    circ.res.means[1,4:5] <- hpd_est_circ(baseline)

    if(length(var.cat) >= 2){

      for(c in seq_along(var.comb.cat[1,])){

        if(!("(Intercept)" %in% colnames(mm$XI))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         rep(0, output$its) + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         output$beta1[, "(Intercept)"] + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         output$beta1[, "(Intercept)"] + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c + length(var.cat)] <- circDiff

        circ.res.cat[c + length(var.cat),1] <- theta_bar(circDiff)
        circ.res.cat[c + length(var.cat),2] <- mode_est_circ(circDiff)
        circ.res.cat[c + length(var.cat),3] <- sd_circ(circDiff)
        circ.res.cat[c + length(var.cat),4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c + 1 + length(var.cat),1] <- theta_bar(dummy)
        circ.res.means[c + 1 + length(var.cat),2] <- mode_est_circ(dummy)
        circ.res.means[c + 1 + length(var.cat),3] <- sd_circ(dummy)
        circ.res.means[c + 1 + length(var.cat),4:5] <- hpd_est_circ(dummy)


    }

    }

  }else{

    for(c in var.cat){

      if(all(!colnames(output$beta1) == c)){

        if(!("(Intercept)" %in% colnames(mm$XI))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + rep(0, output$its))

        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + rep(0, output$its))

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + rep(0, output$its))

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }else if(all(!colnames(output$beta2) == c)){

        if(!("(Intercept)" %in% colnames(mm$XI))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + rep(0, output$its), rep(0, output$its) + output$beta1[, c])

        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + rep(0, output$its), output$beta1[, "(Intercept)"] + output$beta1[, c])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + rep(0, output$its), output$beta1[, "(Intercept)"] + output$beta1[, c])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }else{

        if(!("(Intercept)" %in% colnames(mm$XI))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + output$beta1[, c])

        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }

    }

    circ.res.means[1,1] <- theta_bar(baseline)
    circ.res.means[1,2] <- mode_est_circ(baseline)
    circ.res.means[1,3] <- sd_circ(baseline)
    circ.res.means[1,4:5] <- hpd_est_circ(baseline)

  }



  if(length(var.num) == 0){

    circ.res <- "There are no numeric predictors in the model"
    a.x <- "There are no numeric predictors in the model"
    a.c <- "There are no numeric predictors in the model"
    b.c <- "There are no numeric predictors in the model"
    SAM <- "There are no numeric predictors in the model"
    AS <- "There are no numeric predictors in the model"
    SSDO <- "There are no numeric predictors in the model"

  }else if(!("(Intercept)" %in% colnames(mm$XI)) & !("(Intercept)" %in% colnames(mm$XI))){

    circ.res <- "There is no intercept in the model"
    a.x <- "There is no intercept in the model"
    a.c <- "There is no intercept in the model"
    b.c <- "There is no intercept in the model"
    SAM <- "There is no intercept in the model"
    AS <- "There is no intercept in the model"
    SSDO <- "There is no intercept in the model"

  }else{

    a.x <- matrix(NA, output$its, length(var.num))
    a.c <- matrix(NA, output$its, length(var.num))
    b.c <- matrix(NA, output$its, length(var.num))
    SAM <- matrix(NA, output$its, length(var.num))
    AS <- matrix(NA, output$its, length(var.num))
    SSDO <- matrix(NA, output$its, length(var.num))

    colnames(a.x) <- var.num
    colnames(a.c) <- var.num
    colnames(b.c) <- var.num
    colnames(SAM) <- var.num
    colnames(AS) <- var.num
    colnames(SSDO) <- var.num

    for(v in var.num){

      if(all(!var.num.I == v)){

        if(!("(Intercept)" %in% colnames(mm$XI))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta2[, v])


        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               rep(0, output$its),
                               output$beta2[, v])

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta2[, v])

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.II[v])))

        ASs <- matrix(NA, output$its, length(mm$XII[,1]))

        for(i in seq_along(mm$XII[,1])){
          ASs[,i] <- circest$bc/(1+(circest$bc*(as.numeric(mm$XII[i,v])-circest$ax)))
        }

      }else if(all(!var.num.II == v)){

        if(!("(Intercept)" %in% colnames(mm$XI))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               rep(0, output$its))


        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta1[, v],
                               rep(0, output$its))

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               rep(0, output$its))

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.I[v])))

        ASs <- matrix(NA, output$its, length(mm$XI[,1]))

        for(i in seq_along(mm$XI[,1])){
          ASs[,i] <- circest$bc/(1+(circest$bc*(as.numeric(mm$XI[i,v])-circest$ax)))
        }

      }else{

        if(!("(Intercept)" %in% colnames(mm$XI))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               output$beta2[, v])


        }else if(!("(Intercept)" %in% colnames(mm$XII))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta1[, v],
                               output$beta2[, v])

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               output$beta2[, v])

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.I[v])))

        ASs <- matrix(NA, output$its, length(mm$XI[,1]))

        for(i in seq_along(mm$XI[,1])){
          ASs[,i] <- circest$bc/(1+(circest$bc*(mm$XI[i,v]-circest$ax)))
        }

      }


      a.x[,v] <- circest$ax
      a.c[,v] <- circest$ac
      b.c[,v] <- circest$bc
      SAM[,v] <- SAMs
      AS[,v] <- rowMeans(ASs)
      SSDO[,v] <- circest$SSDO

      circ.res[v,1] <- mean(a.x[,v])
      circ.res[v,2] <- mode_est(a.x[,v])
      circ.res[v,3] <- sd(a.x[,v])
      circ.res[v,4:5] <- hpd_est(a.x[,v])

      circ.res[v,6] <- theta_bar(a.c[,v])
      circ.res[v,7] <- mode_est_circ(a.c[,v])
      circ.res[v,8] <- sd_circ(a.c[,v])
      circ.res[v,9:10] <- hpd_est_circ(a.c[,v])

      circ.res[v,11] <- mean(b.c[,v])
      circ.res[v,12] <- mode_est(b.c[,v])
      circ.res[v,13] <- sd(b.c[,v])
      circ.res[v,14:15] <- hpd_est(b.c[,v])

      circ.res[v,16] <- mean(AS[,v])
      circ.res[v,17] <- mode_est(AS[,v])
      circ.res[v,18] <- sd(AS[,v])
      circ.res[v,19:20] <- hpd_est(AS[,v])

      circ.res[v,21] <- mean(SAM[,v])
      circ.res[v,22] <- mode_est(SAM[,v])
      circ.res[v,23] <- sd(SAM[,v])
      circ.res[v,24:25] <- hpd_est(SAM[,v])

      circ.res[v,26] <- mean(SSDO[,v])
      circ.res[v,27] <- mode_est(SSDO[,v])
      circ.res[v,28] <- sd(SSDO[,v])
      circ.res[v,29:30] <- hpd_est(SSDO[,v])
    }
  }

  model.fit <- DIC_reg(output, mm$XI, mm$XII)


  list(lin.res.I = lin.res.I, lin.res.II = lin.res.II,
       circ.res = circ.res, circ.res.cat = circ.res.cat,
       circ.res.means = circ.res.means,
       a.x = a.x, a.c = a.c, b.c = b.c,
       SAM = SAM, AS = AS, SSDO = SSDO, circ.diff = circ.diff,
       beta1 = output$beta1, beta2 = output$beta2,
       model.fit = model.fit, var.num = var.num, var.cat = var.cat)

}

#' Compute summary and model fit statistics for the circular mixed-effects
#' regression model
#'
#' @param output from a mixed-effects model
#' @param mm output from the function mmme()
#'
#' @keywords internal
#'

summe <- function(output, mm){

  #var.names.I = colnames(mm$mm.I)
  #var.names.II = colnames(mm$mm.II)

  if("(Intercept)" %in% colnames(mm$mm.I)){

    cat <- apply(mm$mm.I, 2, cat_check)
    var.cat.I <- colnames(mm$mm.I)[cat][-1]
    var.num.I <- colnames(mm$mm.I)[!cat]
    grandmeans.I <- colMeans(mm$mm.I)[-1]

    cat.ran <- apply(mm$mm_ran.I, 2, cat_check)
    var.cat.rand.I <- colnames(mm$mm_ran.I)[cat.ran][-1]
    var.num.rand.I <- colnames(mm$mm_ran.I)[!cat.ran]

  }else{

    cat <- apply(mm$mm.I, 2, cat_check)
    var.cat.I <- colnames(mm$mm.I)[cat]
    var.num.I <- colnames(mm$mm.I)[!cat]
    grandmeans.I <- colMeans(mm$mm.I)

    cat.ran <- apply(mm$mm_ran.I, 2, cat_check)
    var.cat.rand.I <- colnames(mm$mm_ran.I)[cat.ran]
    var.num.rand.I <- colnames(mm$mm_ran.I)[!cat.ran]

  }

  if("(Intercept)" %in% colnames(mm$mm.II)){

    cat <- apply(mm$mm.II, 2, cat_check)
    var.cat.II <- colnames(mm$mm.II)[cat][-1]
    var.num.II <- colnames(mm$mm.II)[!cat]
    grandmeans.II <- colMeans(mm$mm.II)[-1]

    cat.ran <- apply(mm$mm_ran.II, 2, cat_check)
    var.cat.rand.II <- colnames(mm$mm_ran.II)[cat.ran][-1]
    var.num.rand.II <- colnames(mm$mm_ran.II)[!cat.ran]

  }else{

    cat <- apply(mm$mm.II, 2, cat_check)
    var.cat.II <- colnames(mm$mm.II)[cat]
    var.num.II <- colnames(mm$mm.II)[!cat]
    grandmeans.II <- colMeans(mm$mm.II)

    cat.ran <- apply(mm$mm_ran.II, 2, cat_check)
    var.cat.rand.II <- colnames(mm$mm_ran.II)[cat.ran]
    var.num.rand.II <- colnames(mm$mm_ran.II)[!cat.ran]

  }

  var.num <- unique(c(var.num.I, var.num.II))
  var.cat <- unique(c(var.cat.I, var.cat.II))

  var.num.rand <- unique(c(var.num.rand.I, var.num.rand.II))
  var.cat.rand <- unique(c(var.cat.rand.I, var.cat.rand.II))

  #Fixed Effects

  lin.res.I <- matrix(NA, ncol(output$beta1), 5)
  rownames(lin.res.I) <- colnames(output$beta1)
  colnames(lin.res.I) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")
  lin.res.II <- matrix(NA, ncol(output$beta2), 5)
  rownames(lin.res.II) <- colnames(output$beta2)
  colnames(lin.res.II) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")

  for(i in colnames(output$beta1)){

    lin.res.I[i,1] <- mean(output$beta1[,i])
    lin.res.I[i,2] <- mode_est(output$beta1[,i])
    lin.res.I[i,3] <- sd(output$beta1[,i])
    lin.res.I[i,4:5] <- hpd_est(output$beta1[,i])

  }

  for(i in colnames(output$beta2)){

    lin.res.II[i,1] <- mean(output$beta2[,i])
    lin.res.II[i,2] <- mode_est(output$beta2[,i])
    lin.res.II[i,3] <- sd(output$beta2[,i])
    lin.res.II[i,4:5] <- hpd_est(output$beta2[,i])

  }

  if(length(var.cat) >= 2){

    var.comb.cat <- combn(var.cat, 2)

    circ.res.means <- matrix(NA, length(var.cat) + ncol(var.comb.cat) + 1, 5)
    names <- sapply(seq_along(var.comb.cat[1,]),
                    function(w){paste(var.comb.cat[,w], sep = "", collapse = "")})
    rownames(circ.res.means) <- c("(Intercept)", var.cat, names)
    colnames(circ.res.means) <- c("mean", "mode", "sd", "LB", "UB")

    circ.diff <- matrix(NA, output$its, length(var.cat) + ncol(var.comb.cat))
    colnames(circ.diff) <- c(var.cat, names)

    circ.res.cat <- matrix(NA, length(var.cat) + ncol(var.comb.cat), 5)
    rownames(circ.res.cat) <- c(var.cat, names)
    colnames(circ.res.cat) <- c("mean", "mode", "sd", "LB", "UB")

  }else{

    var.comb.cat <- 0

    circ.res.means <- matrix(NA, length(var.cat) + 1, 5)
    rownames(circ.res.means) <- c("(Intercept)", var.cat)
    colnames(circ.res.means) <- c("mean", "mode", "sd", "LB", "UB")

    circ.diff <- matrix(NA, output$its, length(var.cat))
    colnames(circ.diff) <- var.cat

    circ.res.cat <- matrix(NA, length(var.cat), 5)
    rownames(circ.res.cat) <- var.cat
    colnames(circ.res.cat) <- c("mean", "mode", "sd", "LB", "UB")

  }



  circ.res <- matrix(NA, length(var.num), 5*6)
  rownames(circ.res) <- var.num
  colnames(circ.res) <- c("mean ax", "mode ax", "sd ax", "LB ax", "UB ax",
                          "mean ac", "mode ac", "sd ac", "LB ac", "UB ac",
                          "mean bc", "mode bc", "sd bc", "LB bc", "UB bc",
                          "mean AS", "mode AS", "sd AS", "LB AS", "UB AS",
                          "mean SAM", "mode SAM", "sd SAM", "LB SAM", "UB SAM",
                          "mean SSDO", "mode SSDO", "sd SSDO", "LB SSSO", "UB SSDO")

  groupmeans.I <- matrix(NA, mm$N, length(var.num.I))
  groupmeans.II <- matrix(NA, mm$N, length(var.num.II))

  colnames(groupmeans.I) <- var.num.I
  colnames(groupmeans.II) <- var.num.II

  if(length(var.cat) == 0){

    circ.res.cat <- "There are no categorical predictors in the model"
    circ.res.means <- "There are no categorical predictors in the model"
    circ.diff <- "There are no categorical predictors in the model"

  }else if(!("(Intercept)" %in% colnames(mm$mm.I)) & !("(Intercept)" %in% colnames(mm$mm.II))){

    circ.res.cat <- "There is no intercept in the model"
    circ.res.means <- "There is no intercept in the model"
    circ.diff <- "There is no intercept in the model"

  }else if(all.equal(var.cat.I, var.cat.II) == TRUE){

    for(c in var.cat){

      if(!("(Intercept)" %in% colnames(mm$mm.I))){

        baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
        dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + output$beta1[, c])

      }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

        baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
        dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

      }else{

        baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
        dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

      }

      diff <- baseline - dummy
      sign <- sign(sin(diff))
      circDiff <- (pi - abs(pi - abs(diff)))*sign

      circ.diff[,c] <- circDiff

      circ.res.cat[c,1] <- theta_bar(circDiff)
      circ.res.cat[c,2] <- mode_est_circ(circDiff)
      circ.res.cat[c,3] <- sd_circ(circDiff)
      circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

      circ.res.means[c,1] <- theta_bar(dummy)
      circ.res.means[c,2] <- mode_est_circ(dummy)
      circ.res.means[c,3] <- sd_circ(dummy)
      circ.res.means[c,4:5] <- hpd_est_circ(dummy)

    }

    circ.res.means[1,1] <- theta_bar(baseline)
    circ.res.means[1,2] <- mode_est_circ(baseline)
    circ.res.means[1,3] <- sd_circ(baseline)
    circ.res.means[1,4:5] <- hpd_est_circ(baseline)

    if(length(var.cat) >= 2){

      for(c in seq_along(var.comb.cat[1,])){

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         rep(0, output$its) + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         output$beta1[, "(Intercept)"] + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, var.comb.cat[1,c]] + output$beta2[, var.comb.cat[2,c]],
                         output$beta1[, "(Intercept)"] + output$beta1[, var.comb.cat[1,c]] + output$beta1[, var.comb.cat[2,c]])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c + length(var.cat)] <- circDiff

        circ.res.cat[c + length(var.cat),1] <- theta_bar(circDiff)
        circ.res.cat[c + length(var.cat),2] <- mode_est_circ(circDiff)
        circ.res.cat[c + length(var.cat),3] <- sd_circ(circDiff)
        circ.res.cat[c + length(var.cat),4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c + 1 + length(var.cat),1] <- theta_bar(dummy)
        circ.res.means[c + 1 + length(var.cat),2] <- mode_est_circ(dummy)
        circ.res.means[c + 1 + length(var.cat),3] <- sd_circ(dummy)
        circ.res.means[c + 1 + length(var.cat),4:5] <- hpd_est_circ(dummy)

      }

    }

  }else{

    for(c in var.cat){

      if(all(!colnames(output$beta1) == c)){

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + rep(0, output$its))

        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + rep(0, output$its))

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + rep(0, output$its))

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }else if(all(!colnames(output$beta2) == c)){

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + rep(0, output$its), rep(0, output$its) + output$beta1[, c])

        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + rep(0, output$its), output$beta1[, "(Intercept)"] + output$beta1[, c])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + rep(0, output$its), output$beta1[, "(Intercept)"] + output$beta1[, c])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }else{

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          baseline <- atan2(output$beta2[, "(Intercept)"], rep(0, output$its))
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], rep(0, output$its) + output$beta1[, c])

        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          baseline <- atan2(rep(0, output$its), output$beta1[, "(Intercept)"])
          dummy <- atan2(rep(0, output$its) + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

        }else{

          baseline <- atan2(output$beta2[, "(Intercept)"], output$beta1[, "(Intercept)"])
          dummy <- atan2(output$beta2[, "(Intercept)"] + output$beta2[, c], output$beta1[, "(Intercept)"] + output$beta1[, c])

        }

        diff <- baseline - dummy
        sign <- sign(sin(diff))
        circDiff <- (pi - abs(pi - abs(diff)))*sign

        circ.diff[,c] <- circDiff

        circ.res.cat[c,1] <- theta_bar(circDiff)
        circ.res.cat[c,2] <- mode_est_circ(circDiff)
        circ.res.cat[c,3] <- sd_circ(circDiff)
        circ.res.cat[c,4:5] <- hpd_est_circ(circDiff)

        circ.res.means[c,1] <- theta_bar(dummy)
        circ.res.means[c,2] <- mode_est_circ(dummy)
        circ.res.means[c,3] <- sd_circ(dummy)
        circ.res.means[c,4:5] <- hpd_est_circ(dummy)

      }

    }

    circ.res.means[1,1] <- theta_bar(baseline)
    circ.res.means[1,2] <- mode_est_circ(baseline)
    circ.res.means[1,3] <- sd_circ(baseline)
    circ.res.means[1,4:5] <- hpd_est_circ(baseline)

  }


  if(length(var.num) == 0){

    circ.res <- "There are no numeric predictors in the model"
    groupmeans.I <- "There are no numeric predictors in the model"
    groupmeans.II <- "There are no numeric predictors in the model"
    grandmeans.I <- "There are no numeric predictors in the model"
    grandmeans.II <- "There are no numeric predictors in the model"
    a.x <- "There are no numeric predictors in the model"
    a.c <- "There are no numeric predictors in the model"
    b.c <- "There are no numeric predictors in the model"
    SAM <- "There are no numeric predictors in the model"
    AS <- "There are no numeric predictors in the model"
    SSDO <- "There are no numeric predictors in the model"

  }else if(!("(Intercept)" %in% colnames(mm$mm.I)) & !("(Intercept)" %in% colnames(mm$mm.II))){

    circ.res <- "There is no intercept in the model"
    groupmeans.I <- "There is no intercept in the model"
    groupmeans.II <- "There is no intercept in the model"
    grandmeans.I <- "There is no intercept in the model"
    grandmeans.II <- "There is no intercept in the model"
    a.x <- "There is no intercept in the model"
    a.c <- "There is no intercept in the model"
    b.c <- "There is no intercept in the model"
    SAM <- "There is no intercept in the model"
    AS <- "There is no intercept in the model"
    SSDO <- "There is no intercept in the model"

  }else{

    a.x <- matrix(NA, output$its, length(var.num))
    a.c <- matrix(NA, output$its, length(var.num))
    b.c <- matrix(NA, output$its, length(var.num))
    SAM <- matrix(NA, output$its, length(var.num))
    AS <- matrix(NA, output$its, length(var.num))
    SSDO <- matrix(NA, output$its, length(var.num))

    colnames(a.x) <- var.num
    colnames(a.c) <- var.num
    colnames(b.c) <- var.num
    colnames(SAM) <- var.num
    colnames(AS) <- var.num
    colnames(SSDO) <- var.num

    for(v in var.num){

      if(all(!var.num.I == v)){

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta2[, v])


        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               rep(0, output$its),
                               output$beta2[, v])

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta2[, v])

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.II[v])))

        for(i in 1:mm$N){
          groupmeans.II[i,v] <- mean(mm$XII[[i]][,v])
        }

        ASs <- matrix(NA, output$its, mm$N)

        for(i in 1:mm$N){
          ASs[,i] <- circest$bc/(1+(circest$bc*(as.numeric(do.call(rbind, mm$XII)[i,v])-circest$ax)))
        }

      }else if(all(!var.num.II == v)){

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               rep(0, output$its))


        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta1[, v],
                               rep(0, output$its))

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               rep(0, output$its))

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.I[v])))

        for(i in 1:mm$N){
          groupmeans.I[i,v] <- mean(mm$XI[[i]][,v])
        }

        ASs <- matrix(NA, output$its, mm$N)

        for(i in 1:mm$N){
          ASs[,i] <- circest$bc/(1+(circest$bc*(as.numeric(do.call(rbind, mm$XI)[i,v])-circest$ax)))
        }

      }else{

        if(!("(Intercept)" %in% colnames(mm$mm.I))){

          circest <- circ_coef(rep(0, output$its),
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               output$beta2[, v])


        }else if(!("(Intercept)" %in% colnames(mm$mm.II))){

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               rep(0, output$its),
                               output$beta1[, v],
                               output$beta2[, v])

        }else{

          circest <- circ_coef(output$beta1[, "(Intercept)"],
                               output$beta2[, "(Intercept)"],
                               output$beta1[, v],
                               output$beta2[, v])

        }

        SAMs <- circest$bc/(1+(circest$bc*-(circest$ax-grandmeans.I[v])))

        for(i in 1:mm$N){
          groupmeans.I[i,v] <- mean(mm$XI[[i]][,v])
        }

        ASs <- matrix(NA, output$its, mm$N)

        for(i in 1:mm$N){
          ASs[,i] <- circest$bc/(1+(circest$bc*(as.numeric(do.call(rbind, mm$XI)[i,v])-circest$ax)))
        }

      }


      a.x[,v] <- circest$ax
      a.c[,v] <- circest$ac
      b.c[,v] <- circest$bc
      SAM[,v] <- SAMs
      AS[,v] <- rowMeans(ASs)
      SSDO[,v] <- circest$SSDO

      circ.res[v,1] <- mean(a.x[,v])
      circ.res[v,2] <- mode_est(a.x[,v])
      circ.res[v,3] <- sd(a.x[,v])
      circ.res[v,4:5] <- hpd_est(a.x[,v])

      circ.res[v,6] <- theta_bar(a.c[,v])
      circ.res[v,7] <- mode_est_circ(a.c[,v])
      circ.res[v,8] <- sd_circ(a.c[,v])
      circ.res[v,9:10] <- hpd_est_circ(a.c[,v])

      circ.res[v,11] <- mean(b.c[,v])
      circ.res[v,12] <- mode_est(b.c[,v])
      circ.res[v,13] <- sd(b.c[,v])
      circ.res[v,14:15] <- hpd_est(b.c[,v])

      circ.res[v,16] <- mean(AS[,v])
      circ.res[v,17] <- mode_est(AS[,v])
      circ.res[v,18] <- sd(AS[,v])
      circ.res[v,19:20] <- hpd_est(AS[,v])

      circ.res[v,21] <- mean(SAM[,v])
      circ.res[v,22] <- mode_est(SAM[,v])
      circ.res[v,23] <- sd(SAM[,v])
      circ.res[v,24:25] <- hpd_est(SAM[,v])

      circ.res[v,26] <- mean(SSDO[,v])
      circ.res[v,27] <- mode_est(SSDO[,v])
      circ.res[v,28] <- sd(SSDO[,v])
      circ.res[v,29:30] <- hpd_est(SSDO[,v])
    }
  }

  #Random Effects

  lin.res.varrand.I <- matrix(NA, output$q1, 5)

  colnames(lin.res.varrand.I) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")
  lin.res.varrand.II <- matrix(NA, output$q2, 5)
  colnames(lin.res.varrand.II) <- c("mean", "mode", "sd", "LB HPD", "UB HPD")

  circ.res.varrand <- matrix(NA, length(var.num.rand)+length(var.cat.rand)+1, 5)
  colnames(circ.res.varrand) <- c("mean", "mode", "sd", "LB", "UB")

  if(output$q1 ==1 & output$q2 ==1){

    rownames(circ.res.varrand) <- c("RI")

  }else{

    rownames(circ.res.varrand) <- c("RI",
                                    paste("RS", c(var.num.rand, var.cat.rand)))

  }

  if(output$q1 == 1){

    rownames(lin.res.varrand.I) <- c("RI")

  }else{
    rownames(lin.res.varrand.I) <- c("RI",
                                     paste("RS",
                                           c(var.num.rand.I, var.cat.rand.I)))
  }

  if(output$q2 == 1){

    rownames(lin.res.varrand.II) <- c("RI")

  }else{
    rownames(lin.res.varrand.II) <- c("RI", paste("RS", c(var.num.rand.II, var.cat.rand.II)))
  }

  for(i in 1:ncol(output$omega1)){

    lin.res.varrand.I[i,1] <- mean(output$omega1[i,i,])
    lin.res.varrand.I[i,2] <- mode_est(output$omega1[i,i,])
    lin.res.varrand.I[i,3] <- sd(output$omega1[i,i,])
    lin.res.varrand.I[i,4:5] <- hpd_est(output$omega1[i,i,])

  }

  for(i in 1:ncol(output$omega2)){

    lin.res.varrand.II[i,1] <- mean(output$omega2[i,i,])
    lin.res.varrand.II[i,2] <- mode_est(output$omega2[i,i,])
    lin.res.varrand.II[i,3] <- sd(output$omega2[i,i,])
    lin.res.varrand.II[i,4:5] <- hpd_est(output$omega2[i,i,])

  }

  circ.varrand.ri <- apply(output$circular.ri, 2, rho_circ)
  circ.res.varrand[1,1] <- mean(1-circ.varrand.ri)
  circ.res.varrand[1,2] <- mode_est(1-circ.varrand.ri)
  circ.res.varrand[1,3] <- sd(1-circ.varrand.ri)
  circ.res.varrand[1,4:5] <- hpd_est(1-circ.varrand.ri)


  circest.rand.num <- array(NA,
                            dim = c(mm$N, output$its, 3),
                            dimnames = list(NULL, NULL, c("ax", "ac", "bc")))

  circest.rand.cat <- array(NA,
                            dim = c(mm$N, output$its, 1),
                            dimnames = list(NULL, NULL, c("circ.diff")))


  if(length(var.num.rand) == 0){

    varrand.num <-  "There are no continuous variables with a random slope"

  }else{

    varrand.num <-  matrix(NA, output$its, length(var.num.rand))
    colnames(varrand.num) <- var.num.rand

    for(v in var.num.rand){

      if(all(!colnames(output$b1) == v)){

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- rep(0, output$its)
          b2 <- output$beta1[, v]+ output$b1[i,v,]

          circest.rand.num[i,,] <- as.matrix(circ_coef(a1, a2, b1, b2))[,1:3]

        }

      }else if(all(!colnames(output$b2) == v)){

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- output$beta1[, v]+ output$b1[i,v,]
          b2 <- rep(0, output$its)

          circest.rand.num[i,,] <- as.matrix(circ_coef(a1, a2, b1, b2))[,1:3]

        }


      }else{

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- output$beta1[, v]+ output$b1[i,v,]
          b2 <- output$beta1[, v]+ output$b1[i,v,]

          circest.rand.num[i,,] <- as.matrix(circ_coef(a1, a2, b1, b2))[,1:3]

        }

      }

      varrand.num[,v] <- apply(circest.rand.num[,,"bc"], 2, var)

      circ.res.varrand[paste("RS", v), 1] <- mean(varrand.num)
      circ.res.varrand[paste("RS", v), 2] <- mode_est(varrand.num)
      circ.res.varrand[paste("RS", v), 3] <- sd(varrand.num)
      circ.res.varrand[paste("RS", v), 4:5] <- hpd_est(varrand.num)


    }

  }

  if(length(var.cat.rand) == 0){

    varrand.cat <- "There are no categorical variables with a random slope"

  }else{

    varrand.cat <- matrix(NA, output$its, length(var.cat.rand))
    colnames(varrand.cat) <- var.cat.rand

    for(c in var.cat.rand){

      if(all(!colnames(output$b1) == c)){

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- rep(0, output$its)
          b2 <- output$beta2[, c]+ output$b2[i,c,]

          baseline <- atan2(a2, a1)
          dummy <- atan2(a2 + b2, a1 + b1)

          diff <- baseline - dummy
          sign <- sign(sin(diff))
          circDiff <- (pi - abs(pi - abs(diff)))*sign

          circest.rand.cat[i,,] <- circDiff

        }

      }else if(all(!colnames(output$b2) == c)){

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- output$beta1[, c]+ output$b1[i,c,]
          b2 <- rep(0, output$its)

          baseline <- atan2(a2, a1)
          dummy <- atan2(a2 + b2, a1 + b1)

          diff <- baseline - dummy
          sign <- sign(sin(diff))
          circDiff <- (pi - abs(pi - abs(diff)))*sign

          circest.rand.cat[i,,] <- circDiff

        }


      }else{

        for(i in 1:mm$N){

          a1 <- output$beta1[,"(Intercept)"] + output$b1[i,1,]
          a2 <- output$beta2[, "(Intercept)"]+ output$b2[i,1,]
          b1 <- output$beta1[, c]+ output$b1[i,c,]
          b2 <- output$beta1[, c]+ output$b1[i,c,]

          baseline <- atan2(a2, a1)
          dummy <- atan2(a2 + b2, a1 + b1)

          diff <- baseline - dummy
          sign <- sign(sin(diff))
          circDiff <- (pi - abs(pi - abs(diff)))*sign

          circest.rand.cat[i,,] <- circDiff

        }

      }

      varrand.cat[,c] <- apply(circest.rand.cat[,,"circ.diff"], 2, rho_circ)

      circ.res.varrand[paste("RS", c), 1] <- mean(1-varrand.cat)
      circ.res.varrand[paste("RS", c), 2] <- mode_est(1-varrand.cat)
      circ.res.varrand[paste("RS", c), 3] <- sd(1-varrand.cat)
      circ.res.varrand[paste("RS", c), 4:5] <- hpd_est(1-varrand.cat)


    }

  }



  #DIC computations
  Q <- output$its
  predictiva.M <- list()
  lppd.inp <- list()
  pWAIC.inp <- list()
  pWAIC2.inp <- list()
  timesum <- matrix(NA, output$its, mm$N)
  timesumM <- rep(NA, mm$N)
  timesum.lppd.inp <- rep(NA, mm$N)
  timesum.pWAIC.inp <- rep(NA, mm$N)
  timesum.pWAIC2.inp <- rep(NA, mm$N)

  if(output$p1 == 1){

    BI <- mode_est(output$beta1)

  }else{

    BI <- apply(output$beta1, 2, mode_est)

  }

  if(output$p2 == 1){

    BII <- mode_est(output$beta2)

  }else{

    BII <- apply(output$beta2, 2, mode_est)
  }

  if(output$q1 > 1  & output$q2 > 1){

    for(i in 1:mm$N){

      predictiva.M[[i]] <- rep(NA, mm$no.Meas[i])
      lppd.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC2.inp[[i]] <- rep(NA, mm$no.Meas[i])

      for(j in 1:mm$no.Meas[i]){

        t.aux     <- mm$theta[[i]][j]
        mu.ij.I   <- c(BI%*%mm$XI[[i]][j,] + apply(output$b1[i,,], 1, mode_est)%*%mm$ZI[[i]][j,])
        mu.ij.II  <- c(BII%*%mm$XII[[i]][j,] + apply(output$b2[i,,], 1, mode_est)%*%mm$ZII[[i]][j,])
        bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
        norm2.M   <- mu.ij.I^2 + mu.ij.II^2
        c <- (1 + ((bb*pnorm(bb))/dnorm(bb)))
        predictiva.M[[i]][j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2.M))*c
        ll.mean.its <- log(mean(output$predictiva[[i]][,j]))
        mean.ll.its <- mean(log(output$predictiva[[i]][,j]))
        pWAIC2.inp[[i]][j] <- var(log(output$predictiva[[i]][,j]))
        lppd.inp[[i]][j] <- ll.mean.its
        pWAIC.inp[[i]][j] <- ll.mean.its - mean.ll.its

      }

      timesum[,i] <- rowSums(log(output$predictiva[[i]]))
      timesumM[i] <- sum(log(predictiva.M[[i]]))
      timesum.lppd.inp[i] <- sum(lppd.inp[[i]][which(is.finite(lppd.inp[[i]]))])
      timesum.pWAIC.inp[i] <- sum(pWAIC.inp[[i]][which(is.finite(pWAIC.inp[[i]]))])
      timesum.pWAIC2.inp[i] <- sum(pWAIC2.inp[[i]][which(is.finite(pWAIC2.inp[[i]]))])

    }

  }else if(output$q1 > 1 & output$q2 == 1){

    for(i in 1:mm$N){

      predictiva.M[[i]] <- rep(NA, mm$no.Meas[i])
      lppd.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC2.inp[[i]] <- rep(NA, mm$no.Meas[i])

      for(j in 1:mm$no.Meas[i]){

        t.aux     <- mm$theta[[i]][j]
        mu.ij.I   <- c(BI%*%mm$XI[[i]][j,] + apply(output$b1[i,,], 1, mode_est)%*%mm$ZI[[i]][j,])
        mu.ij.II  <- c(BII%*%mm$XII[[i]][j,] + apply(t(as.matrix(output$b2[i,,])), 1, mode_est)%*%mm$ZII[[i]][j,])
        bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
        norm2.M   <- mu.ij.I^2 + mu.ij.II^2
        c <- (1 + ((bb*pnorm(bb))/dnorm(bb)))
        predictiva.M[[i]][j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2.M))*c
        ll.mean.its <- log(mean(output$predictiva[[i]][,j]))
        mean.ll.its <- mean(log(output$predictiva[[i]][,j]))
        pWAIC2.inp[[i]][j] <- var(log(output$predictiva[[i]][,j]))
        lppd.inp[[i]][j] <- ll.mean.its
        pWAIC.inp[[i]][j] <- ll.mean.its - mean.ll.its

      }

      timesum[,i] <- rowSums(log(output$predictiva[[i]]))
      timesumM[i] <- sum(log(predictiva.M[[i]]))
      timesum.lppd.inp[i] <- sum(lppd.inp[[i]][which(is.finite(lppd.inp[[i]]))])
      timesum.pWAIC.inp[i] <- sum(pWAIC.inp[[i]][which(is.finite(pWAIC.inp[[i]]))])
      timesum.pWAIC2.inp[i] <- sum(pWAIC2.inp[[i]][which(is.finite(pWAIC2.inp[[i]]))])

    }


  }else if(output$q1 == 1 & output$q2 > 1){

    for(i in 1:mm$N){

      predictiva.M[[i]] <- rep(NA, mm$no.Meas[i])
      lppd.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC2.inp[[i]] <- rep(NA, mm$no.Meas[i])

      for(j in 1:mm$no.Meas[i]){

        t.aux     <- mm$theta[[i]][j]
        mu.ij.I   <- c(BI%*%mm$XI[[i]][j,] + apply(t(as.matrix(output$b1[i,,])), 1, mode_est)%*%mm$ZI[[i]][j,])
        mu.ij.II  <- c(BII%*%mm$XII[[i]][j,] + apply(output$b2[i,,], 1, mode_est)%*%mm$ZII[[i]][j,])
        bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
        norm2.M   <- mu.ij.I^2 + mu.ij.II^2
        c <- (1 + ((bb*pnorm(bb))/dnorm(bb)))
        predictiva.M[[i]][j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2.M))*c
        ll.mean.its <- log(mean(output$predictiva[[i]][,j]))
        mean.ll.its <- mean(log(output$predictiva[[i]][,j]))
        pWAIC2.inp[[i]][j] <- var(log(output$predictiva[[i]][,j]))
        lppd.inp[[i]][j] <- ll.mean.its
        pWAIC.inp[[i]][j] <- ll.mean.its - mean.ll.its

      }

      timesum[,i] <- rowSums(log(output$predictiva[[i]]))
      timesumM[i] <- sum(log(predictiva.M[[i]]))
      timesum.lppd.inp[i] <- sum(lppd.inp[[i]][which(is.finite(lppd.inp[[i]]))])
      timesum.pWAIC.inp[i] <- sum(pWAIC.inp[[i]][which(is.finite(pWAIC.inp[[i]]))])
      timesum.pWAIC2.inp[i] <- sum(pWAIC2.inp[[i]][which(is.finite(pWAIC2.inp[[i]]))])

    }

  }else{

    for(i in 1:mm$N){

      predictiva.M[[i]] <- rep(NA, mm$no.Meas[i])
      lppd.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC.inp[[i]] <- rep(NA, mm$no.Meas[i])
      pWAIC2.inp[[i]] <- rep(NA, mm$no.Meas[i])

      for(j in 1:mm$no.Meas[i]){

        t.aux     <- mm$theta[[i]][j]
        mu.ij.I   <- c(BI%*%mm$XI[[i]][j,] + apply(t(as.matrix(output$b1[i,,])), 1, mode_est)%*%mm$ZI[[i]][j,])
        mu.ij.II  <- c(BII%*%mm$XII[[i]][j,] + apply(t(as.matrix(output$b2[i,,])), 1, mode_est)%*%mm$ZII[[i]][j,])
        bb        <- Dbd(t.aux, mu.ij.I, mu.ij.II)
        norm2.M   <- mu.ij.I^2 + mu.ij.II^2
        c <- (1 + ((bb*pnorm(bb))/dnorm(bb)))
        predictiva.M[[i]][j] <- as.numeric((1/(2*pi))*exp(-0.5*norm2.M))*c
        ll.mean.its <- log(mean(output$predictiva[[i]][,j]))
        mean.ll.its <- mean(log(output$predictiva[[i]][,j]))
        pWAIC2.inp[[i]][j] <- var(log(output$predictiva[[i]][,j]))
        lppd.inp[[i]][j] <- ll.mean.its
        pWAIC.inp[[i]][j] <- ll.mean.its - mean.ll.its

      }

      timesum[,i] <- rowSums(log(output$predictiva[[i]]))
      timesumM[i] <- sum(log(predictiva.M[[i]]))
      timesum.lppd.inp[i] <- sum(lppd.inp[[i]][which(is.finite(lppd.inp[[i]]))])
      timesum.pWAIC.inp[i] <- sum(pWAIC.inp[[i]][which(is.finite(pWAIC.inp[[i]]))])
      timesum.pWAIC2.inp[i] <- sum(pWAIC2.inp[[i]][which(is.finite(pWAIC2.inp[[i]]))])

    }

  }



  personsum <- rowSums(timesum)
  personsum <- personsum[which(is.finite(personsum))]
  personsumM <- sum(timesumM)
  personsumM <- personsumM[which(is.finite(personsumM))]


  lppd <-  sum(timesum.lppd.inp)
  Dhat <- personsumM
  Dbar <- mean(personsum)

  pD <- 2*(Dhat-Dbar)
  pV <- 2*var(personsum)
  pWAIC <- 2*sum(timesum.pWAIC.inp)
  pWAIC2 <- sum(timesum.pWAIC2.inp)

  WAIC <- -2*(lppd - pWAIC)
  WAIC2 <- -2*(lppd - pWAIC2)
  DIC <- -2*(Dhat - pD)
  DICalt <- -2*(Dhat - pV)

  model.fit <- cbind(lppd,
                     Dhat, Dbar, DIC, DICalt, WAIC, WAIC2,
                     pD, pV, pWAIC, pWAIC2)
  colnames(model.fit)<- c("lppd",
                          "Dhat", "Dbar", "DIC", "DICalt", "WAIC", "WAIC2",
                          "pD", "pV", "pWAIC", "pWAIC2")


  list(model.fit = model.fit,
       lin.res.I = lin.res.I,
       lin.res.II = lin.res.II,
       varrand.num = varrand.num,
       varrand.cat = varrand.cat,
       circ.varrand.ri = circ.varrand.ri,
       circ.res = circ.res,
       circ.res.cat = circ.res.cat,
       circ.res.means = circ.res.means,
       a.x = a.x, a.c = a.c, b.c = b.c,
       SAM = SAM, AS = AS, SSDO = SSDO, circ.diff = circ.diff,
       lin.res.varrand.I = lin.res.varrand.I,
       lin.res.varrand.II = lin.res.varrand.II,
       circ.res.varrand = circ.res.varrand)
}
