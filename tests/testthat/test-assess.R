# assess_one() decides what each posterior is compared against. Two ways it can
# be wrong without erroring: mapping a parameter to the wrong truth, and
# crediting a parameter that was never actually estimated.

post_row <- function(parameter, mean, sd = 0.05, lo = mean - 0.1, hi = mean + 0.1) {
  tibble(parameter = parameter, mean = mean, sd = sd, lo = lo, hi = hi,
         max_rhat = 1, min_ess_bulk = 1000, n_divergent = 0L,
         n_max_treedepth = 0L, min_ebfmi = 1)
}

base_true <- function(...) {
  tibble(treatment_effect_mean_normalised = -0.333,
         treatment_effect_sd_normalised   = 0.067,
         time_trend_mean_normalised       = -0.089,
         time_trend_sd_normalised         = 0.044,
         baseline_mean = 0.45, ...)
}

test_that("parameters are scored against the scale the fit actually used", {
  p <- post_row("treatment_effect_mean", -0.33)
  tr <- base_true(treatment_effect_mean_raw = -0.15)
  norm <- assess_one(p, tr, list(normalise = TRUE))
  raw  <- assess_one(p, tr, list(normalise = FALSE))
  expect_equal(norm$true_value, -0.333)
  expect_equal(raw$true_value,  -0.15)
})

test_that("a pinned parameter is dropped, not credited with perfect recovery", {
  # kappa fixed, or mu_gamma pinned at zero, gives sd == 0. Scoring that would
  # report flawless recovery of a constant.
  p <- bind_rows(post_row("treatment_effect_mean", -0.33),
                 post_row("kappa", 0.5, sd = 0))
  out <- assess_one(p, base_true(kappa = 1), list(normalise = TRUE))
  expect_false("kappa" %in% out$parameter)
  expect_true("treatment_effect_mean" %in% out$parameter)
})

test_that("a genuinely sampled kappa is scored", {
  p <- bind_rows(post_row("treatment_effect_mean", -0.33),
                 post_row("kappa", 0.9, sd = 0.3))
  out <- assess_one(p, base_true(kappa = 1), list(normalise = TRUE))
  expect_true("kappa" %in% out$parameter)
  expect_equal(out$true_value[out$parameter == "kappa"], 1)
})

test_that("tau_gamma targets sqrt(mu^2 + tau^2) when mu_gamma is pinned at zero", {
  # With the population mean pinned, a real one-sided imbalance has nowhere to
  # go but the spread. Scoring tau alone would mark correct behaviour a failure
  # -- exactly what the zero-mean and two-population arms do in the transport
  # scenarios.
  tr <- base_true(baseline_difference_mean_normalised = 0.18,
                  baseline_difference_sd_normalised   = 0.049)
  p  <- post_row("baseline_difference_sd", 0.187, sd = 0.02)
  pinned <- assess_one(p, tr, list(normalise = TRUE, mu_gamma = "zero"))
  free   <- assess_one(p, tr, list(normalise = TRUE, mu_gamma = "estimated"))
  expect_equal(pinned$true_value, sqrt(0.18^2 + 0.049^2))
  expect_equal(free$true_value,   0.049)
})

test_that("the tau_gamma adjustment is inert when there is no mean to absorb", {
  tr <- base_true(baseline_difference_mean_normalised = 0,
                  baseline_difference_sd_normalised   = 0.05)
  out <- assess_one(post_row("baseline_difference_sd", 0.05), tr,
                    list(normalise = TRUE, mu_gamma = "zero"))
  expect_equal(out$true_value, 0.05)
})

test_that("parameters with no truth available are dropped silently", {
  # build_true_params() withholds gamma truths when the populations disagree;
  # assess_one() must not invent one.
  p <- bind_rows(post_row("treatment_effect_mean", -0.33),
                 post_row("baseline_difference_sd", 0.05))
  out <- assess_one(p, base_true(), list(normalise = TRUE))
  expect_false("baseline_difference_sd" %in% out$parameter)
})

test_that("coverage and bias are computed from the interval, not the point", {
  p  <- post_row("treatment_effect_mean", -0.30, lo = -0.40, hi = -0.20)
  out <- assess_one(p, base_true(), list(normalise = TRUE))
  expect_true(out$covers)                       # -0.333 lies inside
  expect_equal(out$bias, -0.30 - (-0.333))
  p2 <- post_row("treatment_effect_mean", -0.10, lo = -0.15, hi = -0.05)
  out2 <- assess_one(p2, base_true(), list(normalise = TRUE))
  expect_false(out2$covers)
})
