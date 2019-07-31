# bpnreg 1.0.0
This is the first version of the package.

## bpnreg 1.0.0.9000
This is a developement version.

The function bpnme can now also handle the standard R data.frame instead of only Data Frame Tbl (from dplyr, tbl_df(data))

The function coef.bpnr now correctly returns radians instead of degrees when option units == "radians" is chosen.

The documentation now states that the circular outcome needs to be measured in radians.

The references in the documentation have been updated.

The functions mmr and mmme now give out a warning message when the circular outcome is out of range.
