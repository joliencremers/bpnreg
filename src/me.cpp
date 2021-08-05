/*
# ------------------------------------------------------------------------------
# me.cpp
#
# Runs an MCMC sampler for a circular mixed-effects model.
#
# ------------------------------------------------------------------------------
*/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>

#include "utils.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

const double pi = boost::math::constants::pi<double>();

//' Sample fixed effect coefficients
//'
//' @param Omega current covariance matrix
//' @param Y current outcome vector (Y.I=cos(theta)*R, Y.II=sin(theta)*R)
//' @param X design matrix model parameters (differs per person)
//' @param Z design matrix for random effects (differs per person)
//' @param p dimension X (number of columns/variables+indicator variables)
//' @param A prior variance of fixed effect coefficients
//' @param N sample size at second level
//'
//' @keywords internal
//'
//[[Rcpp::export]]

arma::mat betaBlock(arma::mat Omega, Rcpp::List R, Rcpp::List theta, Rcpp::List X, Rcpp::List Z, int p, arma::mat A, int N){

  //compute inverse of precision matrix (so this is the covariance matrix)
  arma::mat invOmega = arma::inv(Omega);

  //compute the covariance and precision (inverse) matrix of the outcome vector

  arma::mat sXtVX = arma::mat(p,p, fill::zeros);
  arma::mat sXtVY_tmp = arma::mat(p, N);
  arma::mat sXtVY;

  for (int i = 0; i < N; ++i){

    arma::mat Z_tmp = Z[i];
    arma::mat X_tmp = X[i];
    arma::mat theta_tmp = theta[i];
    arma::mat R_temp = R[i];
    arma::mat Y = theta_tmp%R_temp;

    int l = Z_tmp.n_rows;
    arma::mat Vi    = Z_tmp*invOmega*Z_tmp.t() + arma::mat(l,l,fill::eye);
    arma::mat invVi = arma::inv(Vi);
    arma::mat aux  = X_tmp.t()*invVi*X_tmp;
    sXtVX = sXtVX + aux;
    sXtVY_tmp.col(i) = X_tmp.t()*invVi*Y;

  }

  //Compute variance of the distribution for the coefficients vector
  arma::mat beta_var = arma::inv(A + sXtVX);

  //Compute the mean of the distribution for the coefficients vector.

  sXtVY = arma::sum(sXtVY_tmp,1);

  if(p > 1){

    arma::colvec beta_mu = beta_var*sXtVY;
    //Sample the coefficients vector
    arma::mat beta = arma::mvnrnd(beta_mu, beta_var);
    return beta;

  }else{

    arma::mat beta_mu = beta_var*sXtVY;
    //Sample the coefficients vector
    arma::mat beta = arma::randn(1)*sqrt(beta_var)+beta_mu;
    return beta;

  }

}


//' Sample subject specific random effects
//'
//' @param Omega current covariance matrix
//' @param beta current fixed effect coefficients vector
//' @param Y current outcome vector (Y.I=cos(theta)*R, Y.II=sin(theta)*R)
//' @param X design matrix model parameters (differs per person)
//' @param q dimension Z(number of random effects)
//' @param Z design matrix for random effects (differs per person)
//' @param ZtZ transpose(Z)*Z
//' @param N sample size at second level
//'
//' @keywords internal
//'
//[[Rcpp::export]]

arma::mat b_samp(arma::mat Omega, arma::mat beta,
                 Rcpp::List R, Rcpp::List theta, Rcpp::List X,
                 int q, Rcpp::List Z, Rcpp::List ZtZ, int N){

  arma::mat b(N, q);

  for (int i = 0; i < N; ++i){

    arma::mat Z_tmp = Z[i];
    arma::mat ZtZ_tmp = ZtZ[i];
    arma::mat X_tmp = X[i];
    arma::mat theta_tmp = theta[i];
    arma::mat R_temp = R[i];
    arma::mat Y = theta_tmp%R_temp;

    //Compute variance of the distribution for the subject specific random effects
    arma::mat invD = arma::inv(ZtZ_tmp + Omega);
    //Compute mean of the distribution for the subject specific random effects
    arma::mat etilde = Y - (X_tmp*beta);
    arma::mat Zetilde = Z_tmp.t()*etilde;
    arma::colvec b_mu = invD*Zetilde;
    //Sample subject specific random effects vectors
    b.row(i) = arma::mvnrnd(b_mu, invD).t();

  }

  return(b);

}

//' Sample precision matrix
//'
//' @param b subject specific random effects vectors
//' @param B prior sum of squares matrix, scale parameter Wishart distribution
//'   (of b=random effect if prior close to 0-->no random effect)
//' @param q dimension Z(number of random effects)
//' @param v prior df=dimension Z
//' @param N sample size at second level
//'
//' @keywords internal
//'
//[[Rcpp::export]]

arma::mat omega_samp(arma::mat b, arma::mat B, int v, int q, int N){

  arma::mat bb = b.t()*b;
  arma::mat Bbb = arma::inv(B + bb);
  arma::mat Omega = arma::wishrnd(Bbb, v + N);

  return(Omega);

}

//' Compute the Likelihood of the PN distribution (mixed effects)
//'
//' @param theta_cos A List with the cosine of the circular dependent variable.
//' @param theta_sin A List with the sine of the circular dependent variable.
//' @param X1 A list of fixed effect model matrices for component I.
//' @param X2 A list of fixed effect model matrices for component II.
//' @param Z1 A list of random effect model matrices for component I.
//' @param Z2 A list of random effect model matrices for component II.
//' @param beta1 estimated fixed effect coefficients of the first component
//' @param beta2 estimated fixed effect coefficients of the second component
//' @param b1 estimated random effect coefficients of the first component
//' @param b2 estimated random effect coefficients of the second component
//' @param N sample size at second level
//' @param pred An empty list for likelihood computation.
//' @param iteration iteration number at which likelihood is computed

//'
// [[Rcpp::export]]

Rcpp::List lik_me(Rcpp::List theta_cos, Rcpp::List theta_sin,
                  Rcpp::List X1, Rcpp::List X2,
                  Rcpp::List Z1, Rcpp::List Z2,
                  arma::mat beta1, arma::mat beta2, arma::mat b1, arma::mat b2,
                  int N, Rcpp::List pred, int iteration){

  for (int i = 0; i < N; ++i){

    arma::mat mtmp = pred[i];

    arma::mat Z1_tmp = Z1[i];
    arma::mat X1_tmp = X1[i];
    arma::mat theta_cos_tmp = theta_cos[i];

    arma::mat Z2_tmp = Z2[i];
    arma::mat X2_tmp = X2[i];
    arma::mat theta_sin_tmp = theta_sin[i];

    arma::rowvec mub1 = beta1.t()*X1_tmp.t() + b1.row(i)*Z1_tmp.t();
    arma::rowvec mub2 = beta2.t()*X2_tmp.t() + b2.row(i)*Z2_tmp.t();
    arma::rowvec Dbd = theta_cos_tmp.t()%mub1 + theta_sin_tmp.t()%mub2;
    arma::rowvec norm2 = arma::pow(mub1, 2) + arma::pow(mub2, 2);

    int n = mub1.n_cols;
    arma::vec c(n, fill::zeros);
    arma::vec L(n, fill::zeros);

    for (int j = 0; j < n; ++j){

      double Dbd_tmp = Dbd.row(0)(j);

      c(j) = 1 + ((Dbd_tmp*R::pnorm(Dbd_tmp, 0, 1, TRUE, FALSE))/R::dnorm(Dbd_tmp,0, 1, FALSE));
      L(j) = (1/(2*pi))*exp(-0.5*norm2.row(0)(j))*c(j);

    }

    mtmp.row(iteration) = L.t();
    pred[i] = mtmp;

  }

  return(pred);

}

//' A Gibbs sampler for a projected normal mixed-effects model
//'
//' @param theta_cos A List with the cosine of the circular dependent variable.
//' @param theta_sin A List with the sine of the circular dependent variable.
//' @param X1 A list of fixed effect model matrices for component I.
//' @param X2 A list of fixed effect model matrices for component II.
//' @param Z1 A list of random effect model matrices for component I.
//' @param Z2 A list of random effect model matrices for component II.
//' @param ZtZ1 A list of transformed random effect model matrices for component I.
//' @param ZtZ2 A list of transformed random effect model matrices for component II.
//' @param R A list of starting values for R.
//' @param pred An empty list for likelihood computation.
//' @param its An integer specifying the number of iterations
//' @param lag An integer specifying the amount of lag.
//' @param burn An integer specifying the number of burn-in iterations.
//' @param N An integer specifying the number of burn-in iterations.
//'
// [[Rcpp::export]]
Rcpp::List pnme(List theta_cos, List theta_sin,
         List X1, List X2, List Z1, List Z2, List ZtZ1, List ZtZ2,
         List R, Rcpp::List pred,
         int its, int lag, int burn,
         int N) {

  int burn_new = burn * lag;
  int tm = burn_new + (its * lag);

  arma::mat X1_temp = X1[0];
  arma::mat X2_temp = X2[0];
  arma::mat Z1_temp = Z1[0];
  arma::mat Z2_temp = Z2[0];

  int p1 = X1_temp.n_cols;
  int p2 = X2_temp.n_cols;
  int q1 = Z1_temp.n_cols;
  int q2 = Z2_temp.n_cols;

  //Specify priors
  int v1 = q1;
  int v2 = q2;

  arma::mat B1 = arma::eye(v1, v1)*0.0001;
  arma::mat B2 = arma::eye(v2, v2)*0.0001;
  arma::mat A1 = arma::eye(p1, p1)*0.0001;
  arma::mat A2 = arma::eye(p2, p2)*0.0001;

  //Create results objects

  arma::mat beta1(p1,its);
  arma::mat beta2(p2,its);
  arma::cube b1(N,q1,its);
  arma::cube b2(N,q2,its);
  arma::cube omega1(q1,q1,its);
  arma::cube omega2(q2,q2,its);

  arma::mat beta1_tmp(p1, 1, fill::zeros);
  arma::mat beta2_tmp(p2, 1, fill::zeros);
  arma::mat b1_tmp(N,q1, fill::zeros);
  arma::mat b2_tmp(N,q1, fill::zeros);

  //Set starting values

  arma::mat omega1_tmp(q1,q1, fill::eye);
  arma::mat omega2_tmp(q2,q2, fill::eye);

  Rcpp::List R_tmp = R;

  for (int it = 0; it < tm ; ++ it){

    beta1_tmp = betaBlock(omega1_tmp, theta_cos, R_tmp, X1, Z1, p1, A1, N);
    beta2_tmp = betaBlock(omega2_tmp, theta_sin, R_tmp, X2, Z2, p2, A2, N);
    b1_tmp = b_samp(omega1_tmp, beta1_tmp, theta_cos, R_tmp, X1, q1, Z1, ZtZ1, N);
    b2_tmp = b_samp(omega2_tmp, beta2_tmp, theta_sin, R_tmp, X2, q2, Z2, ZtZ2, N);
    omega1_tmp = omega_samp(b1_tmp, B1, v1, q1, N);
    omega2_tmp = omega_samp(b2_tmp, B2, v2, q2, N);

    //Sample latent lengths
    for (int i = 0; i < N; ++i){

      arma::vec R_temp = R_tmp[i];

      arma::mat Z1_tmp = Z1[i];
      arma::mat X1_tmp = X1[i];
      arma::mat theta_cos_tmp = theta_cos[i];

      arma::mat Z2_tmp = Z2[i];
      arma::mat X2_tmp = X2[i];
      arma::mat theta_sin_tmp = theta_sin[i];

      arma::rowvec mub1 = beta1_tmp.t()*X1_tmp.t() + b1_tmp.row(i)*Z1_tmp.t();
      arma::rowvec mub2 = beta2_tmp.t()*X2_tmp.t() + b2_tmp.row(i)*Z2_tmp.t();
      arma::rowvec b = theta_cos_tmp.t()%mub1 + theta_sin_tmp.t()%mub2;

      for (int j = 0; j < R_temp.n_elem; ++j){

        arma::vec y = as<arma::vec>(runif(1,0,exp(-.5*pow(R_temp(j)-b.row(0)(j), 2))));

        double r1 = b.row(0)(j) + max(-b(j), -sqrt(-2*log(y(0))));
        double r2 = b.row(0)(j) + sqrt(-2*log(y(0)));
        R_temp(j)  = sqrt(((pow(r2,2)-pow(r1,2)) * arma::randu()) + pow(r1,2));

      }

    R_tmp[i] = R_temp;

    }




    if ((it + 1 - burn_new > 0) & ((it + 1 -burn_new) % lag == 0)){

      int ii = (it + 1 - burn_new) / lag;
      Rcout << "Iteration:" << ii << "\n";

      beta1.col(ii-1) = beta1_tmp;
      beta2.col(ii-1) = beta2_tmp;
      b1.slice(ii-1) = b1_tmp;
      b2.slice(ii-1) = b2_tmp;
      omega1.slice(ii-1) = omega1_tmp;
      omega2.slice(ii-1) = omega2_tmp;

      pred = lik_me(theta_cos, theta_sin, X1, X2, Z1, Z2,
                    beta1_tmp, beta2_tmp, b1_tmp, b2_tmp,
                    N, pred, ii-1);

    }

  }

  return  Rcpp::List::create(Rcpp::Named("beta1") = beta1,
                             Rcpp::Named("beta2") = beta2,
                             Rcpp::Named("b1") = b1,
                             Rcpp::Named("b2") = b2,
                             Rcpp::Named("omega1") = omega1,
                             Rcpp::Named("omega2") = omega2,
                             Rcpp::Named("its") = its,
                             Rcpp::Named("n.lag") = lag,
                             Rcpp::Named("burn-in") = burn_new,
                             Rcpp::Named("p1") = p1,
                             Rcpp::Named("p2") = p2,
                             Rcpp::Named("q1") = q1,
                             Rcpp::Named("q2") = q2,
                             Rcpp::Named("predictiva") = pred);

}
