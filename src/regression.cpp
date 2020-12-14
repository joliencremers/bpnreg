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

arma::vec lik_reg(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat beta1, arma::mat beta2, int n){

  arma::mat mub1 = beta1*X1.t();
  arma::mat mub2 = beta2*X2.t();
  arma::mat Dbd = cos(theta)%mub1.t() + sin(theta)%mub2.t();
  arma::mat norm2 = arma::pow(mub1, 2) + arma::pow(mub2, 2);
  arma::vec L(n);
  arma::vec c(n);

  for (int jjj=0; jjj < n; ++jjj){
    c(jjj) = 1 + ((Dbd.row(jjj)(0)*R::pnorm(Dbd.row(jjj)(0), 0, 1, TRUE, FALSE))/R::dnorm(Dbd.row(jjj)(0),0, 1, FALSE));
    L(jjj) = (1/(2*pi))*exp(-0.5*norm2.col(jjj)(0))*c(jjj);
  }

  return L;

}

//' Compute Model Fit Measures Regression Model
//'
//' @param Output output from the circular regression function Regression()
//' @param X1 model matrix for the first component
//' @param X2 model matrix for the second component
//'
//[[Rcpp::export]]

Rcpp::List DIC_reg(Rcpp::List output, arma::mat X1, arma::mat X2){

  arma::vec theta = Rcpp::as<arma::vec>(output["theta"]);
  arma::mat beta1 = Rcpp::as<arma::mat>(output["beta1"]);
  arma::mat beta2 = Rcpp::as<arma::mat>(output["beta2"]);
  arma::mat Likelihood = Rcpp::as<arma::mat>(output["Likelihood"]);

  double n = theta.n_elem;
  double p1 = X1.n_cols;
  double p2 = X2.n_cols;

  arma::mat u = arma::join_rows(cos(theta),sin(theta));

  arma::mat m_beta1 = arma::mean(beta1, 0);
  arma::mat m_beta2 = arma::mean(beta2, 0);

  arma::vec Likelihood_m = lik_reg(X1, X2, theta, m_beta1, m_beta2, n);

  arma::mat personsum = arma::sum(log(Likelihood), 1);
  double D_hat = arma::sum(log(Likelihood_m));
  arma::vec D_bar = arma::mean(personsum,0);

  double pD = 2*(D_hat - D_bar(0));
  arma::vec pV = 2*var(personsum);
  double pV_alt = 2*(p1 + p2);

  double DIC = -2*D_hat + 2*pD;
  double DIC_alt = -2*D_hat + 2*pV(0);
  double DIC_alt2 = -2*D_hat + 2*pV_alt;

  arma::mat lppd = arma::sum(log(arma::mean(Likelihood, 0)),1);

  arma::mat pWAIC = 2*arma::sum(log(arma::mean(Likelihood, 0)) - arma::mean(log(Likelihood), 0),1);
  arma::mat pWAIC2 = arma::sum(var(log(Likelihood), 0),1);

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

arma::mat slice_rcpp(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat beta1, arma::mat beta2, int n, arma::mat r){

  arma::mat mub1 = beta1*X1.t();
  arma::mat mub2 = beta2*X2.t();
  arma::mat Dbd = cos(theta)%mub1.t() + sin(theta)%mub2.t();
  for (int jjj=0; jjj < n; ++jjj){

    arma::mat y = as<arma::vec>(runif(1,0,1)) % exp(-.5*pow((r.row(jjj)-Dbd.row(jjj)),2));
    arma::mat u = as<arma::vec>(runif(1,0,1));

    arma::mat r1 = Dbd.row(jjj) + max(-Dbd.row(jjj), -sqrt(-2*log(y)));
    arma::mat r2 = Dbd.row(jjj) + sqrt(-2*log(y));
    r.row(jjj)  = sqrt(((pow(r2,2)-pow(r1,2)) % u) + pow(r1,2));
  }
  return r;
}


//' A Gibbs sampler for a projected normal regression model
//'
//' @param theta A vector with the circular dependent variable.
//' @param X1r A model matrix for component I.
//' @param X2r A model matrix for component II.
//' @param its An integer specifying the number of iterations
//' @param lag An integer specifying the amount of lag.
//' @param burn An integer specifying the number of burn-in iterations.
//'
// [[Rcpp::export]]

Rcpp::List pnr(arma::vec theta,
               NumericMatrix X1r,
               NumericMatrix X2r,
               int its, int lag, int burn) {

  arma::mat X1 = Rcpp::as<arma::mat>(X1r);
  arma::mat X2 = Rcpp::as<arma::mat>(X2r);

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

  int bb = burn*lag;
  int kk = bb + (its*lag);

  arma::mat r = arma::ones<arma::mat>(n,1);
  arma::mat beta1_tmp(kk, p1);
  arma::mat beta2_tmp(kk, p2);
  arma::mat predictiva_tmp(kk, n);

  arma::mat beta1(its, p1);
  arma::mat beta2(its, p2);
  arma::mat predictiva(its, n);

  arma::mat Y = r%datose.each_col();

  //Gibbs iterations
  for(int iii=0; iii < kk; ++iii){

    arma::mat XtY1 = X1.t()*Y.col(0);
    arma::mat XtY2 = X2.t()*Y.col(1);
    arma::mat mstar1 = sigma1*(v1mu1 + XtY1);
    arma::mat mstar2 = sigma2*(v2mu2 + XtY2);

    //Sample coefficients
    beta1_tmp.row(iii) = mvrnorm_arma_eigen(1, mstar1.col(0), sigma1).t();
    beta2_tmp.row(iii) = mvrnorm_arma_eigen(1, mstar2.col(0), sigma2).t();

    //Sample R
    r = slice_rcpp(X1, X2, theta, beta1_tmp.row(iii), beta2_tmp.row(iii), n, r);

    //Compute Y
    Y = r%datose.each_col();

    predictiva_tmp.row(iii) = lik_reg(X1, X2, theta, beta1_tmp.row(iii), beta2_tmp.row(iii), n).t();

  }

  for(int i=0; i < its; ++i){

    int index = (i*lag) + bb;

    beta1.row(i) = beta1_tmp.row(index);
    beta2.row(i) = beta2_tmp.row(index);
    predictiva.row(i) = predictiva_tmp.row(index);

  }

  return  Rcpp::List::create(Rcpp::Named("beta1") = beta1,
                             Rcpp::Named("beta2") = beta2,
                             Rcpp::Named("Likelihood") = predictiva,
                             Rcpp::Named("its") = its,
                             Rcpp::Named("lag") = lag,
                             Rcpp::Named("burn-in") = burn,
                             Rcpp::Named("p1") = p1,
                             Rcpp::Named("p2") = p2,
                             Rcpp::Named("theta") = theta);

}


