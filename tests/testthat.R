# Fast, deterministic tests for the simulation scaffolding.
#
# These deliberately fit NOTHING. In a simulation study the expensive part is
# the MCMC, but the DANGEROUS part is the deterministic scaffolding around it:
# a wrong truth value, a mis-mapped parameter or a DGP that does not generate
# what the scenario asked for does not error -- it produces a plausible,
# confident, wrong number. Those are what this suite covers, in seconds.
#
# Run with:  testthat::test_dir("tests/testthat")
library(testthat)
test_dir("tests/testthat")
