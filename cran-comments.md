## Submission

This is a patch to fix valgrind and macOS errors from version 2.0.0, submitted
to CRAN on 02-03-2021.

In this version I have:

* Removed MASS from imports to fix note in the 2.0.0 version

* Changed c++ code to solve errors from CRAN valgrind checks in the 2.0.0 version.

* Changed arma::inv_sympd to arma::inv in c++ code to solve errors on macOS in the 2.0.0 version.

## Test environments:

*Local Windows Install, R 4.0.1

*CRAN win-builder (devel and release)

*Ubuntu 16.04.6 (devel and release)

*Debian Linux with valgrind (release) 

*macOS 10.13.6 High Sierra on rhub (release)

## R CMD check results: 

There were no ERRORs or WARNINGs on the development version and current release of R.

## Downstream dependencies:

There are currently no downstream dependencies for this package.
