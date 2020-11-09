data <- bpnreg::Maps

pred.I = Error.rad ~ Maze + Trial.type + (1|Subject)
pred.II = pred.I
its = 1000
burn = 1
n.lag = 1
seed = 101


if (!is.null(seed)){set.seed(seed)}

mm <- bpnreg:::mmme(pred.I, data, pred.II)


X1r <- mm$XI
X2r <- mm$XII
Z1r <- mm$ZI
Z2r <- mm$ZII
theta <- mm$theta


Rcpp::sourceCpp('C:/Users/bsj777/Onderzoek/R-packages/bpnreg/src/me.cpp')

pnme(theta, X1r, X2r, Z1r, Z2r, its, n.lag, burn)
