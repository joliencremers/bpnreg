## Resubmission This is a resubmission. In this version I have:

* Removed capitalization of Normal in the description file except in the title 
which has to be title case

* Included references for the method in the 'Description' field of the
DESCRIPTION file.

## Test environments: 

*Local Windows Install, R 3.4.1

*CRAN win-builder (devel and release)

*Ubuntu 14.04 (on travis.ci), R 3.4.2

## R CMD check results There were no ERRORs or WARNINGs

There were 2 NOTES:

* checking CRAN incoming feasibility ... NOTE Maintainer: 'Jolien Cremers 
<joliencremers@gmail.com>'

Note used as reminder for CRAN maintainers, safe to ignore.

* Checking compiled code ... NOTE File 'bpnreg/libs/x64/bpnreg.dll': Found no 
calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol search.

See 'Writing portable packages' in the 'Writing R Extensions' manual.

This is a known false positive on Windows. It did not occur on the CRAN 
win-builder or on Ubuntu.

## Downstream dependencies 

There are currently no downstream dependencies for this package.
