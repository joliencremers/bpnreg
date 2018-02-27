/*
# ------------------------------------------------------------------------------
# regression.cpp
#
# Includes functions to compute a circular mean and a mean resultant length
# Runs an MCMC sampler for a circular regression model.
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

NumericVector circ_coef_rcpp(double a1, double a2, double b1, double b2) {

  double ax = -(a1*b1 + a2*b2)/(pow(b1,2) + pow(b2,2));
  double ac = atan2(a2 + b2*ax, a1 + b1*ax);
  double bc = (tan(atan2(a2, a1))-ac/-ax);
  double SDO = sqrt(pow(a1 + b1, 2) + pow(a2 + b2, 2));
  double SSDO = R::sign(sin(ac - atan2(b2, b1)))*SDO;

  Rcpp::NumericVector out = NumericVector::create(ax, ac, bc, SDO, SSDO);

  out.names() = CharacterVector::create("ax", "ac", "bc", "SDO", "SSDO");

  return out;
}


//' Compute the Likelihood of the PN distribution
//'
//' @param X1 the model matrix of the first component
//' @param X2 the model matrix of the second component
//' @param theta a circular outcome value
//' @param b1 estimated linear coefficients of the first component
//' @param b2 estimated linear coefficients of the second component
//' @param n sample size
//'
//[[Rcpp::export]]

arma::vec lik(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat b1, arma::mat b2, int n){

  arma::mat mub1 = b1*X1.t();
  arma::mat mub2 = b2*X2.t();
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

Rcpp::List DIC_reg(Rcpp::List Output, arma::mat X1, arma::mat X2){

  arma::vec theta = Rcpp::as<arma::vec>(Output["theta"]);
  arma::mat B1 = Rcpp::as<arma::mat>(Output["B1"]);
  arma::mat B2 = Rcpp::as<arma::mat>(Output["B2"]);
  arma::mat Likelihood = Rcpp::as<arma::mat>(Output["Likelihood"]);

  double n = theta.n_elem;
  double p1 = X1.n_cols;
  double p2 = X2.n_cols;

  arma::mat u = arma::join_rows(cos(theta),sin(theta));

  arma::mat m_B1 = arma::mean(B1, 0);
  arma::mat m_B2 = arma::mean(B2, 0);

  arma::vec Likelihood_m = lik(X1, X2, theta, m_B1, m_B2, n);

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
//' @param b1 A matrix containing the coefficients of component I for the current iteration.
//' @param b2 A matrix containing the coefficients of component II for the current iteration.
//' @param n An integer indicating the sample size of the data.
//' @param r A matrix with the estimates of r of the previous iteration.
//'
//[[Rcpp::export]]

arma::mat slice_rcpp(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat b1, arma::mat b2, int n, arma::mat r){

  arma::mat mub1 = b1*X1.t();
  arma::mat mub2 = b2*X2.t();
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
  arma::mat B1(kk, p1);
  arma::mat B2(kk, p2);
  arma::mat bc(kk, p1-1);
  arma::mat SSDO(kk, p1-1);
  arma::mat ac(kk, p1-1);
  arma::mat ax(kk, p1-1);
  arma::mat SAM(kk, p1-1);
  arma::mat AS(kk, p1-1);
  arma::mat Predictiva(kk, n);

  arma::mat B1s(its, p1);
  arma::mat B2s(its, p2);
  arma::mat bcs(its, p1-1);
  arma::mat SSDOs(its, p1-1);
  arma::mat acs(its, p1-1);
  arma::mat axs(its, p1-1);
  arma::mat SAMs(its, p1-1);
  arma::mat ASs(its, p1-1);
  arma::mat Predictivas(its, n);

  arma::mat Y = r%datose.each_col();

  //Gibbs iterations
  for(int iii=0; iii < kk; ++iii){

    arma::mat XtY1 = X1.t()*Y.col(0);
    arma::mat XtY2 = X2.t()*Y.col(1);
    arma::mat mstar1 = sigma1*(v1mu1 + XtY1);
    arma::mat mstar2 = sigma2*(v2mu2 + XtY2);

    //Sample coefficients
    arma::mat b1 = mvrnorm_arma_eigen(1, mstar1.col(0), sigma1).t();
    arma::mat b2 = mvrnorm_arma_eigen(1, mstar2.col(0), sigma2).t();

    //Sample R
    r = slice_rcpp(X1, X2, theta, b1, b2, n, r);

    //Compute Y
    Y = r%datose.each_col();

    //Fill posterior coefficient matrices
    B1.row(iii) = b1;
    B2.row(iii) = b2;

    Predictiva.row(iii) = lik(X1, X2, theta, b1, b2, n).t();

  }

  for(int i=0; i < its; ++i){

    int index = (i*lag) + bb;

    B1s.row(i) = B1.row(index);
    B2s.row(i) = B2.row(index);
    Predictivas.row(i) = Predictiva.row(index);
  }

  return  Rcpp::List::create(Rcpp::Named("B1") = B1s,
                             Rcpp::Named("B2") = B2s,
                             Rcpp::Named("Likelihood") = Predictivas,
                             Rcpp::Named("its") = its,
                             Rcpp::Named("lag") = lag,
                             Rcpp::Named("burn-in") = burn,
                             Rcpp::Named("p1") = p1,
                             Rcpp::Named("p2") = p2,
                             Rcpp::Named("theta") = theta);

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

NumericVector hmodeciC(NumericVector x, double cip) {

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

NumericVector hmodeci(NumericVector x, double cip) {

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
