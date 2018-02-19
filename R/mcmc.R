#' Sample fixed effect coefficients in a Random Intercept model
#'
#' @param Omega current covariance matrix
#' @param Y current outcome vector (Y.I=cos(theta)*R, Y.II=sin(theta)*R)
#' @param X design matrix model parameters (differs per person)
#' @param Z design matrix for random effects (differs per person)
#' @param p dimension X (number of columns/variables+indicator variables)
#' @param A prior variance of fixed effect coefficients
#' @param N sample size at second level

betaBlock.fRI <- function(Omega, Y, X, Z, p, A, N){

  #compute inverse of precision matrix (so this is the covariance matrix)
  invOmega <- solve(Omega)

  #compute the covariance and precision (inverse) matrix of the outcome vector
  Vi    <- list(NA)
  invVi <- list(NA)

  for (i in 1:N){
    Vi[[i]]    <- (Z[[i]] %*% invOmega %*% t(Z[[i]])) + diag(length(Z[[i]]))
    invVi[[i]] <- chol2inv(chol(Vi[[i]]))
  }

  #Compute variance of the distribution for the coefficients vector
  sXtVX <- diag(p) * 0.0
  aux   <- lapply( (1:N), function(w){t(X[[w]]) %*% invVi[[w]] %*% X[[w]]} )

  for (i in 1:N){sXtVX <- sXtVX + aux[[i]]}

  Var.beta <- chol2inv(chol(A + sXtVX ))

  #Compute the mean of the distribution for the coefficients vector.
  sXtVY <- vapply( (1:N),
                   function(w){t(X[[w]]) %*% invVi[[w]] %*% Y[[w]]},
                   FUN.VALUE = rep(1,p))

  if(p > 1){
    sXtVY <- rowSums(sXtVY)
  }else{
    sXtVY <- sum(sXtVY)
  }

  betaF <- c(Var.beta %*% sXtVY)


  #Sample the coefficients vector.
  beta.aux <- MASS::mvrnorm(1, betaF, Var.beta)

  drop(beta.aux)
}

#' Sample fixed effect coefficients in a Random Slope model
#'
#' @inheritParams betaBlock.fRI

betaBlock.fRS <- function(Omega, Y, X, Z, p, A, N){

  #compute inverse of precision matrix (so this is the covariance matrix)
  invOmega <- solve(Omega)

  #compute the covariance and precision (inverse) matrix of the outcome vector
  Vi    <- list(NA)
  invVi <- list(NA)

  for (i in 1:N){
    Vi[[i]]    <- (Z[[i]] %*% invOmega %*% t(Z[[i]])) + diag(length(Z[[i]][,1]))
    invVi[[i]] <- chol2inv(chol(Vi[[i]]))
  }

  #Compute variance of the distribution for the coefficients vector
  sXtVX <- diag(p) * 0.0
  aux   <- lapply( (1:N), function(w){t(X[[w]]) %*% invVi[[w]] %*% X[[w]]} )

  for (i in 1:N){sXtVX <- sXtVX + aux[[i]]}

  Var.beta <- chol2inv(chol(A + sXtVX ))

  #Compute the mean of the distribution for the coefficients vector.
  sXtVY <-vapply( (1:N),
                  function(w){t(X[[w]]) %*% invVi[[w]] %*% Y[[w]]},
                  FUN.VALUE = rep(1,p))

  if (p > 1){
    sXtVY <- rowSums(sXtVY)
  }else{
    sXtVY <- sum(sXtVY)
  }

  betaF <- c(Var.beta %*% sXtVY)

  #Sample the coefficients vector.
  beta.aux<-MASS::mvrnorm(1, betaF, Var.beta)

  drop(beta.aux)
}

#' Sample subject specific random effects
#'
#' @param Omega current covariance matrix
#' @param beta current fixed effect coefficients vector
#' @param Y current outcome vector (Y.I=cos(theta)*R, Y.II=sin(theta)*R)
#' @param X design matrix model parameters (differs per person)
#' @param q dimension Z(number of random effects)
#' @param Z design matrix for random effects (differs per person)
#' @param ZtZ transpose(Z)*Z
#' @param N sample size at second level

b.f<-function(Omega, beta, Y, X, q, Z, ZtZ, N){

  bF   <- array(NA, c(N, q))
  D    <- array(NA, c(q, q, N))
  invD <- array(NA, c(q, q, N))

  #Compute variance of the distribution for the subject specific random effects
  for (i in 1:N){
    D[,,i]    <- (ZtZ[[i]]) + Omega
    invD[,,i] <- solve(D[,,i])
  }

  #Compute mean of the distribution for the subject specific random effects

  for (i in 1:N){
    etilde <- Y[[i]] - (X[[i]] %*% beta)
    bF[i,] <- t(invD[,,i] %*% t(Z[[i]]) %*% etilde)
  }

  #Sample subject specific random effects vectors
  b.aux <- t(as.matrix(vapply( (1:N),
                               function(w){MASS::mvrnorm(1, bF[w,], invD[,,w])},
                               FUN.VALUE = rep(1,q))))

  drop(b.aux)
}

#' Sample precision matrix
#'
#' @param b subject specific random effects vectors
#' @param B prior sum of squares matrix, scale parameter Wishart distribution
#'   (of b=random effect if prior close to 0-->no random effect)
#' @param q dimension Z(number of random effects)
#' @param v prior df=dimension Z
#' @param N sample size at second level

Omega.f<-function(b, B, v, q, N){

  Omatrix <- matrix(0, q, q)
  bb.sum  <- t(b) %*% b
  Bbb     <- B + bb.sum
  vc      <- solve(Bbb)
  Omatrix <- rWishart(1, v + N, vc)

  drop(Omatrix)
}

#' Sample R (latent lengths)
#'
#' @param t current outcome value
#' @param mu1 current predicted linear mean of the first component
#' @param mu2 current predicted linear mean of the second component
#' @param r current value for r

slice_r_me<-function(t, mu1, mu2, r){
  b    <- Dbd(t, mu1, mu2)
  y    <- runif(1, 0, exp(-0.5 * (r-b)^2) )
  u    <- runif(1, 0, 1)
  r1   <- b + max(c(-b, -sqrt(-2 * log(y))))
  r2   <- b + sqrt(-2 * log(y))
  rnew <- sqrt( ( (r2^2 - r1^2) * u) + r1^2)

  drop(rnew)
}

#' Compute utmu
#'
#' @param t current outcome value
#' @param mu1 current predicted linear mean of the first component
#' @param mu2 current predicted linear mean of the second component

Dbd <- function(t, mu1, mu2){ cos(t)*mu1 + sin(t)*mu2 }
