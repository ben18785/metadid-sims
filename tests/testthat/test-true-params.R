# build_true_params() computes the values every bias and coverage number in the
# study is scored against. A bug here does not error -- it silently redefines
# what "correct" means, and the whole report reads as a finding.
#
# This is not hypothetical: R/simulate.R records that scoring the per-study
# model against E[theta]/E[baseline] rather than E[theta_i/b_i] is what
# produced the spurious "G8 bias" (ben18785/metadid#39).

test_that("normalised truths follow the documented delta-method expansion", {
  dgp <- mk_dgp(true_effect = -0.15, sigma_effect = 0.03,
                true_trend  = -0.04, sigma_trend  = 0.02,
                baseline_mean = 0.45, baseline_sd = 0.05)
  tp <- build_true_params(dgp, list(), normalised = TRUE)

  bm <- 0.45; cv_b2 <- (0.05 / 0.45)^2

  expect_equal(tp$treatment_effect_mean_normalised, (-0.15 / bm) * (1 + cv_b2))
  expect_equal(tp$time_trend_mean_normalised,       (-0.04 / bm) * (1 + cv_b2))
  expect_equal(tp$treatment_effect_sd_normalised,
               sqrt(0.03^2 + (-0.15)^2 * cv_b2) / bm)
  expect_equal(tp$time_trend_sd_normalised,
               sqrt(0.02^2 + (-0.04)^2 * cv_b2) / bm)
})

test_that("with no between-study baseline variation the expansion is exact", {
  # CV_b = 0 removes the Jensen term, so the normalised truth must collapse to
  # the naive raw/bm. If it does not, the correction is being applied twice.
  dgp <- mk_dgp(true_effect = -0.15, sigma_effect = 0.03, baseline_sd = 0)
  tp  <- build_true_params(dgp, list(), normalised = TRUE)
  expect_equal(tp$treatment_effect_mean_normalised, -0.15 / 0.45)
  expect_equal(tp$treatment_effect_sd_normalised,    0.03 / 0.45)
})

test_that("raw truths are untouched by the normalisation machinery", {
  dgp <- mk_dgp(true_effect = -0.15, sigma_effect = 0.03, baseline_sd = 0.05)
  tp  <- build_true_params(dgp, list(), normalised = FALSE)
  expect_equal(tp$treatment_effect_mean_raw, -0.15)
  expect_equal(tp$treatment_effect_sd_raw,    0.03)
  expect_equal(tp$baseline_mean,              0.45)
})

# --- gamma truths ----------------------------------------------------------
# The model carries ONE non-randomised gamma population. When the DiD and RCT
# populations differ, no true value exists and emitting one would score the
# model against a number we invented.

test_that("gamma truths are emitted when the non-randomised populations agree", {
  dgp <- mk_dgp(n_did = 10L, n_rct = 10L,
                did_gamma_mean = 0.08, did_gamma_sd = 0.02,
                rct_gamma_mean = 0.08, rct_gamma_sd = 0.02)
  tp <- build_true_params(dgp, list(), normalised = TRUE)
  expect_true("baseline_difference_mean_normalised" %in% names(tp))
  expect_equal(tp$baseline_difference_mean_raw, 0.08)
  expect_equal(tp$baseline_difference_sd_raw,   0.02)
})

test_that("gamma truths are withheld when the populations disagree", {
  dgp <- mk_dgp(n_did = 10L, n_rct = 10L,
                did_gamma_mean = 0.08, did_gamma_sd = 0.02,
                rct_gamma_mean = 0.01, rct_gamma_sd = 0.01)
  tp <- build_true_params(dgp, list(), normalised = TRUE)
  expect_false("baseline_difference_mean_normalised" %in% names(tp))
  expect_false("baseline_difference_sd_normalised"   %in% names(tp))
})

test_that("a design with no non-randomised studies does not veto the truth", {
  # All RCTs randomised, so only the DiD population is non-randomised; its
  # values are the truth and disagreement with the (absent) RCT population is
  # irrelevant.
  dgp <- mk_dgp(n_did = 10L, n_rct = 10L, n_randomised_rct = 10L,
                did_gamma_mean = 0.08, did_gamma_sd = 0.02,
                rct_gamma_mean = 0.99, rct_gamma_sd = 0.99)
  tp <- build_true_params(dgp, list(), normalised = TRUE)
  expect_equal(tp$baseline_difference_mean_raw, 0.08)
})

test_that("a sign-flipped population has mean zero and the pooled magnitude", {
  # Direction randomised study by study: the mean vanishes, the magnitude does
  # not. Getting this wrong would make the direction scenarios unscoreable.
  dgp <- mk_dgp(n_did = 10L, n_rct = 0L,
                did_gamma_mean = 0.08, did_gamma_sd = 0.02,
                did_gamma_sign_flip = TRUE)
  tp <- build_true_params(dgp, list(), normalised = TRUE)
  expect_equal(tp$baseline_difference_mean_raw, 0)
  expect_equal(tp$baseline_difference_sd_raw, sqrt(0.08^2 + 0.02^2))
})

test_that("kappa truth appears only when some study is randomised", {
  expect_false("kappa" %in% names(
    build_true_params(mk_dgp(n_did = 10L, kappa_true = 1), list(), TRUE)))
  tp <- build_true_params(
    mk_dgp(n_did = 10L, n_randomised_did = 10L, kappa_true = 1), list(), TRUE)
  expect_equal(tp$kappa, 1)
})

test_that("design offsets propagate to both scales", {
  dgp <- mk_dgp(true_effect = -0.15, baseline_sd = 0)
  tp  <- build_true_params(dgp, list(delta_rct = 0.02, delta_pp = -0.03), TRUE)
  expect_equal(tp$treatment_effect_mean_rct_raw, -0.13)
  expect_equal(tp$treatment_effect_mean_pp_raw,  -0.18)
  expect_equal(tp$treatment_effect_mean_rct_normalised, -0.13 / 0.45)
})
