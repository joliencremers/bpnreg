/*
# ------------------------------------------------------------------------------
# regression.cpp
#
# Includes functions to run an MCMC sampler for a circular regression model.
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

//' Compute the Likelihood of the PN distribution (regression)
//'
//' @param X1 the model matrix of the first component
//' @param X2 the model matrix of the second component
//' @param theta a circular outcome value
//' @param beta1 estimated linear coefficients of the first component
//' @param beta2 estimated linear coefficients of the second component
//' @param n sample size
//'
//[[Rcpp::export]]

arma::vec lik_reg(arma::mat X1, arma::mat X2, arma::vec theta, arma::vec beta1, arma::vec beta2, int n){

  arma::vec mub1 = X1*beta1;
  arma::vec mub2 = X2*beta2;
  arma::vec Dbd = cos(theta)%mub1 + sin(theta)%mub2;

  arma::vec norm2 = arma::pow(mub1, 2) + arma::pow(mub2, 2);
  arma::vec L(n, fill::zeros);
  arma::vec c(n, fill::zeros);

  for (int jjj=0; jjj < n; ++jjj){
    c(jjj) = 1 + ((Dbd(jjj)*R::pnorm(Dbd(jjj), 0, 1, TRUE, FALSE))/R::dnorm(Dbd(jjj),0, 1, FALSE));
    L(jjj) = (1/(2*pi))*exp(-0.5*norm2(jjj))*c(jjj);
  }

  return L;

}

//' Compute Model Fit Measures Regression Model
//'
//' @param theta circular outcome values
//' @param beta1 regression coefficients for the second component for each mcmc iteration from pnr function
//' @param beta2 regression coefficients for the second component for each mcmc iteration from pnr function
//' @param Likelihood likelihood values for each individual and mcmc itertion from pnr function
//' @param X1 model matrix for the first component
//' @param X2 model matrix for the second component
//'
//[[Rcpp::export]]

Rcpp::List DIC_reg(arma::mat theta, arma::mat beta1, arma::mat beta2, arma::mat Likelihood, arma::mat X1, arma::mat X2){

  double n = theta.n_rows;
  double p1 = X1.n_cols;
  double p2 = X2.n_cols;

  arma::mat u = arma::join_rows(cos(theta),sin(theta));

  arma::rowvec m_beta1 = arma::mean(beta1, 0);
  arma::rowvec m_beta2 = arma::mean(beta2, 0);

  arma::vec Likelihood_m = lik_reg(X1, X2, theta, m_beta1.t(), m_beta2.t(), n);

  arma::vec personsum = arma::sum(trunc_log(Likelihood), 1);
  double D_hat = arma::sum(log(Likelihood_m));
  double D_bar = arma::mean(personsum);

  double pD = 2*(D_hat - D_bar);
  double pV = 2*var(personsum);
  double pV_alt = 2*(p1 + p2);

  double DIC = -2*D_hat + 2*pD;
  double DIC_alt = -2*D_hat + 2*pV;
  double DIC_alt2 = -2*D_hat + 2*pV_alt;

  arma::mat lppd = arma::sum(log(arma::mean(Likelihood, 0)),1);

  arma::mat pWAIC = 2*arma::sum(log(arma::mean(Likelihood, 0)) - arma::mean(trunc_log(Likelihood), 0),1);
  arma::mat pWAIC2 = arma::sum(var(trunc_log(Likelihood), 0),1);

  double WAIC = -2*(lppd.col(0)(0) - pWAIC.col(0)(0));
  double WAIC2 = -2*(lppd.col(0)(0) - pWAIC2.col(0)(0));

  return  Rcpp::List::create(Rcpp::Named("pD") = pD,
                             Rcpp::Named("pV") = pV,
                             Rcpp::Named("pV_alt") = pV_alt,
                             Rcpp::Named("pWAIC") = pWAIC,
                             Rcpp::Named("pWAIC2") = pWAIC2,
                             Rcpp::Named("lppd") = lppd,
                             Rcpp::Named("D_hat") = D_hat,
                             Rcpp::Named("D_bar") = D_bar,
                             Rcpp::Named("DIC") = DIC,
                             Rcpp::Named("DIC_alt") = DIC_alt,
                             Rcpp::Named("DIC_alt2") = DIC_alt2,
                             Rcpp::Named("WAIC") = WAIC,
                             Rcpp::Named("WAIC2") = WAIC2);

}

//' A slice sampler for the latent lengths r
//'
//' @param X1 A model matrix for component I.
//' @param X2 A model matrix for component II.
//' @param theta A vector with the circular dependent variable.
//' @param beta1 A matrix containing the coefficients of component I for the current iteration.
//' @param beta2 A matrix containing the coefficients of component II for the current iteration.
//' @param n An integer indicating the sample size of the data.
//' @param r A matrix with the estimates of r of the previous iteration.
//'
//[[Rcpp::export]]

arma::vec slice_rcpp(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat beta1, arma::mat beta2, int n, arma::vec r){

  arma::vec mub1 = X1*beta1;
  arma::vec mub2 = X2*beta2;
  arma::vec Dbd = cos(theta)%mub1 + sin(theta)%mub2;

  for (int jjj=0; jjj < n; ++jjj){

    double y = R::runif(0,1)*exp(-.5*pow(r(jjj)-Dbd(jjj),2));
    double u = R::runif(0,1);

    double r1 = Dbd(jjj) + max(-Dbd(jjj), -sqrt(-2*log(y)));
    double r2 = Dbd(jjj) + sqrt(-2*log(y));
    r(jjj)  = sqrt(((pow(r2,2)-pow(r1,2))*u) + pow(r1,2));
  }
  return r;
}


//' A Gibbs sampler for a projected normal regression model
//'
//' @param theta A vector with the circular dependent variable.
//' @param X1 A model matrix for component I.
//' @param X2 A model matrix for component II.
//' @param its An integer specifying the number of iterations
//' @param lag An integer specifying the amount of lag.
//' @param burn An integer specifying the number of burn-in iterations.
//'
// [[Rcpp::export]]

Rcpp::List pnr(arma::vec theta,
               arma::mat X1,
               arma::mat X2,
               int its, int lag, int burn) {

  double n = theta.n_elem;
  double p1 = X1.n_cols;
  double p2 = X2.n_cols;

  arma::mat datose(n, 2);
  datose.col(0) = cos(theta);
  datose.col(1) = sin(theta);

  //Prior specification regression parameters
  arma::vec mu1 = arma::ones<arma::vec>(p1)*0;
  arma::vec mu2 = arma::ones<arma::vec>(p2)*0;
  arma::mat v1 = arma::eye<arma::mat>(p1,p1)*0.0001;
  arma::mat v2 = arma::eye<arma::mat>(p2,p2)*0.0001;

  //Posterior specification regression parameters
  arma::mat XtX1 = X1.t()*X1;
  arma::mat XtX2 = X2.t()*X2;
  arma::mat vstar1 = v1 + XtX1;
  arma::mat vstar2 = v2 + XtX2;
  arma::mat sigma1 = arma::inv(vstar1);
  arma::mat sigma2 = arma::inv(vstar2);
  arma::mat v1mu1 = v1*mu1;
  arma::mat v2mu2 = v2*mu2;

  //Initialize matrices for results

  int burn_new = burn*lag;
  int tm = burn_new + (its*lag);

  arma::mat r = arma::ones<arma::mat>(n,1);

  arma::mat beta1_tmp;
  arma::mat beta2_tmp;
  arma::mat predictiva_tmp;

  arma::mat beta1(its, p1);
  arma::mat beta2(its, p2);
  arma::mat predictiva(its, n);

  arma::mat Y = r%datose.each_col();

  //Gibbs iterations
  for(int it=0; it < tm; ++it){

    arma::mat XtY1 = X1.t()*Y.col(0);
    arma::mat XtY2 = X2.t()*Y.col(1);
    arma::mat mstar1 = sigma1*(v1mu1 + XtY1);
    arma::mat mstar2 = sigma2*(v2mu2 + XtY2);

    //Sample coefficients
    beta1_tmp = mvrnorm_arma_eigen(1, mstar1.col(0), sigma1).t();
    beta2_tmp = mvrnorm_arma_eigen(1, mstar2.col(0), sigma2).t();

    //Sample R
    r = slice_rcpp(X1, X2, theta, beta1_tmp.t(), beta2_tmp.t(), n, r);

    //Compute Y
    Y = r%datose.each_col();

    predictiva_tmp = lik_reg(X1, X2, theta, beta1_tmp.t(), beta2_tmp.t(), n).t();

    if ((it + 1 - burn_new > 0) & ((it + 1 -burn_new) % lag == 0)){

      int ii = (it + 1 - burn_new) / lag;
      Rcout << "Iteration:" << ii << "\n";

      beta1.row(ii-1) = beta1_tmp;
      beta2.row(ii-1) = beta2_tmp;
      predictiva.row(ii-1) = predictiva_tmp;

    }

  }

  return  Rcpp::List::create(Rcpp::Named("beta1") = beta1,
                             Rcpp::Named("beta2") = beta2,
                             Rcpp::Named("Likelihood") = predictiva,
                             Rcpp::Named("its") = its,
                             Rcpp::Named("n.lag") = lag,
                             Rcpp::Named("burn-in") = burn,
                             Rcpp::Named("p1") = p1,
                             Rcpp::Named("p2") = p2,
                             Rcpp::Named("theta") = theta);

}


