/*
# ------------------------------------------------------------------------------
# utils.cpp
#
# Includes utility functions for PN regression and mixed-effects models
#
# ------------------------------------------------------------------------------
 */


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

const double pi = boost::math::constants::pi<double>();


//' Compute a mean resultant length
//'
//' @param theta a circular variable in radians.
//'
// [[Rcpp::export]]

Rcpp::List rho(arma::vec theta){

  double n = theta.n_elem;

  double S = arma::sum(arma::sin(theta));
  double C = arma::sum(arma::cos(theta));

  double rho = sqrt(pow(S, 2) + pow(C, 2))/n;

  return  Rcpp::List::create(Rcpp::Named("rho") = rho,
                             Rcpp::Named("C") = C,
                             Rcpp::Named("S") = S);

}


//' Compute a mean direction
//'
//' @inheritParams rho
//'
// [[Rcpp::export]]

double theta_bar(arma::vec theta){

  Rcpp::List R_bar = rho(theta);

  double C = R_bar["C"];
  double S = R_bar["S"];
  double rho = R_bar["rho"];

  double mean_dir = atan2(S/rho, C/rho);

  return mean_dir;

}


//' Compute Eigenvalues
//'
//' @param X A matrix.
//'
// [[Rcpp::export]]

arma::vec eigen_val(arma::mat X) {

  arma::cx_vec eigval;
  arma::eig_gen(eigval, X);

  return arma::conv_to<arma::vec>::from(eigval);
}

//' Compute Eigenvectors
//'
//' @inheritParams eigen_val
//'
// [[Rcpp::export]]

arma::mat eigen_vec(arma::mat X) {

  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, X);

  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

//' Sample from a multivariate normal distribution
//'
//' @param sigma A variance-covariance matrix.
//' @param mu A mean vector.
//' @param n An integer indicating the number of samples to take.
//'
// [[Rcpp::export]]

arma::mat mvrnorm_arma_eigen(int n, arma::vec mu, arma::mat sigma) {

  int ncols = sigma.n_cols;
  arma::mat Y(n, ncols);

  for(int iii=0; iii < ncols; ++iii){
    Y.col(iii) = as<arma::vec>(rnorm(n));
  }

  return arma::repmat(mu, 1, n) + (eigen_vec(sigma) * arma::diagmat(sqrt(max(eigen_val(sigma), arma::zeros(ncols)))) * Y.t());
}

//' Compute circular coefficients
//'
//' @param a1 intercept estimate of component I.
//' @param a2 intercept estimate of component I.
//' @param b1 slope estimate of component I.
//' @param b2 slope estimate of component I.
//'
// [[Rcpp::export]]

Rcpp::NumericVector circ_coef_rcpp(double a1, double a2, double b1, double b2) {

  double ax = -(a1*b1 + a2*b2)/(pow(b1,2) + pow(b2,2));
  double ac = atan2(a2 + b2*ax, a1 + b1*ax);
  double bc = (tan(atan2(a2, a1))-ac/-ax);
  double SDO = sqrt(pow(a1 + b1, 2) + pow(a2 + b2, 2));
  double SSDO = R::sign(sin(ac - atan2(b2, b1)))*SDO;

  Rcpp::NumericVector out = NumericVector::create(ax, ac, bc, SDO, SSDO);

  out.names() = CharacterVector::create("ax", "ac", "bc", "SDO", "SSDO");

  return out;
}

//' Estimate the mode by finding the highest posterior density interval
//'
//' @param x a  sample from which to estimate the interval
//' @param cip bandwidth for the algorithm, ranging from 0 to 1
//'
//' @return a scalar containing the estimate of the mode
//'
// [[Rcpp::export]]

double hmodeC(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln, M;

  n = x.size();
  NumericVector sx = clone(x);
  NumericVector sx2 = clone(x)+(2*pi);
  std::vector<double> SX;
  SX.reserve( x.size() + x.size() ); // preallocate memory
  SX.insert( SX.end(), sx.begin(), sx.end() );
  SX.insert( SX.end(), sx2.begin(), sx2.end() );
  std::sort(SX.begin(), SX.end());
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = SX[cil]-SX[0];

  for (int i=0; i < (n); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (SX[i+cil]-SX[i])) {
      ln = (SX[i+cil]-SX[i]);
      chiv = i;
    }
  }

  M = (fmod(SX[chiv+cil],(2*pi))+SX[chiv])/2;

  return M;
}



//' Find the highest density interval of a circular variable
//'
//' @inheritParams hmodeC
//'
//' @return a vector of length 2 containing the lower and upper bound of the interval
//'
// [[Rcpp::export]]

Rcpp::NumericVector hmodeciC(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln;

  n = x.size();
  NumericVector sx = clone(x);
  NumericVector sx2 = clone(x)+(2*pi);
  std::vector<double> SX;
  SX.reserve( x.size() + x.size() ); // preallocate memory
  SX.insert( SX.end(), sx.begin(), sx.end() );
  SX.insert( SX.end(), sx2.begin(), sx2.end() );
  std::sort(SX.begin(), SX.end());
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Length of the currently smallest interval.
  ln = SX[cil]-SX[0];

  for (int i=0; i < (n); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (SX[i+cil]-SX[i])) {
      ln = (SX[i+cil]-SX[i]);
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = SX[chiv];
  M[1] = fmod(SX[chiv+cil],(2*pi));

  return M;
}


//' Estimate the mode by finding the highest posterior density interval
//'
//' @inheritParams hmodeC
//'
//' @return a scalar containing the estimate of the mode
//'
// [[Rcpp::export]]

double hmode(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln, M;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  M = (sx[chiv+cil]+sx[chiv])/2;

  return M;
}


//' Find the highest density interval.
//'
//' @inheritParams hmodeC
//' @inheritParams hmodeC
//'
//' @return a vector of length 2 containing the lower and upper bound of the interval.
//'
// [[Rcpp::export]]

Rcpp::NumericVector hmodeci(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Length of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = sx[chiv];
  M[1] = sx[chiv+cil];

  return M;
}
