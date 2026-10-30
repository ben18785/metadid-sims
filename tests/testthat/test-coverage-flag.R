# coverage_flag() replaces a fixed threshold, which could not work across the
# replication counts the scenarios run at. These pin the thresholds it implies,
# because a silent change to them would relabel defects as noise or the reverse.

test_that("the flagging thresholds adapt to the replication count", {
  # Flag at or below: 0.667 (n=15), 0.720 (n=25), 0.778 (n=45), 0.800 (n=50).
  expect_true(coverage_flag(10 / 15, 15))
  expect_false(coverage_flag(11 / 15, 15))     # the A12 / I7 case: one replicate
  expect_true(coverage_flag(18 / 25, 25))      # the B6 zero_mean case
  expect_false(coverage_flag(19 / 25, 25))
  expect_true(coverage_flag(35 / 45, 45))
  expect_false(coverage_flag(36 / 45, 45))
  expect_true(coverage_flag(40 / 50, 50))
  expect_false(coverage_flag(41 / 50, 50))
})

test_that("a calibrated scenario is rarely flagged, a broken one usually is", {
  # The false-alarm rate is the significance level by construction, and well
  # under the 5.6% that a flat 0.80 threshold produced at n = 15.
  set.seed(1)
  for (n in c(15L, 25L, 50L)) {
    cov_ok <- rbinom(4000, n, 0.90) / n
    expect_lt(mean(coverage_flag(cov_ok, rep(n, 4000))), 0.03)
  }
  # And it still catches a scenario that really is off, where there are enough
  # replicates to tell.
  cov_bad <- rbinom(4000, 50L, 0.75) / 50
  expect_gt(mean(coverage_flag(cov_bad, rep(50L, 4000))), 0.7)
})

test_that("missing or zero-replication rows are not flagged", {
  expect_equal(coverage_flag(c(NA, 0.5, 0.5), c(25, NA, 0)), c(FALSE, FALSE, FALSE))
})
