# Does the simulated data actually carry the parameters it was asked for?
#
# This is the most fundamental check in the repo and the easiest to skip. Every
# bias and coverage number assumes the DGP put theta, beta and the baselines
# where the scenario said. If it did not, the study measures something other
# than what it reports, and the posterior gives no hint -- the fits converge
# happily on the wrong data.
#
# Structural parameters are checked from the RAW simulated means, independently
# of build_true_params(), so a matching bug in both cannot cancel out.
#
# Heterogeneity is switched off (sigma_effect = sigma_trend = baseline_sd = 0)
# so the target is exact rather than a realised draw, and n is large so
# sampling error is small. Tolerances come from the analytic SE, not from trial
# and error: with S studies, n per arm, within-SD s and correlation rho, the
# double-difference has SE = sqrt(4 * s^2 * (1 - rho) / n / S). At S = 400,
# n = 200, s = 0.12 and rho = 0.5 that is ~0.0006, so TOL is ~10 SE.
#
# The generous margin is deliberate. These tests are deterministic -- fixed
# seeds -- so they cannot flake as written, but a refactor that shifts how much
# randomness is consumed re-rolls every draw. A seed sitting at 3 SE would pass
# today and fail for no real reason later. Ten SE of headroom still catches a
# theta error of ~0.01 (7% of the effect), verified by perturbing the target.
#
# The simulator itself is unbiased here: across ten seeds the mean deviation of
# the double-difference from theta is +0.00001 (SE 0.00037), with and without a
# non-zero gamma.
S_STUDIES <- 400L
TOL <- 0.006

exact_dgp <- function(...) {
  mk_dgp(sigma_effect = 0, sigma_trend = 0, baseline_sd = 0,
         within_sd = 0.12, rho = 0.5,
         n_control = 200L, n_treatment = 200L, ...)
}

# --- simulate_imbalance_grid (categories X) --------------------------------

grid_summary <- function(..., seed = 202L) {
  cfg <- scenario("t", dgp = c(list(type = "bespoke",
                                    bespoke_fn = "simulate_imbalance_grid"),
                               as.list(exact_dgp(...))))
  simulate_scenario("t", seed, cfg)$data$summary_data
}

test_that("the DiD double-difference recovers theta", {
  d <- grid_summary(n_did = S_STUDIES, n_rct = 0L,
                    true_effect = -0.15, true_trend = -0.04,
                    did_gamma_mean = 0.08, did_gamma_sd = 0.02)
  dd <- (d$mean_post_treatment - d$mean_pre_treatment) -
        (d$mean_post_control   - d$mean_pre_control)
  # The double difference must return theta REGARDLESS of gamma -- that is the
  # whole point of the design, and it is what makes gamma a nuisance rather
  # than a confounder for DiD.
  expect_equal(mean(dd), -0.15, tolerance = TOL)
})

test_that("the control arm's change recovers the time trend", {
  d <- grid_summary(n_did = S_STUDIES, n_rct = 0L,
                    true_effect = -0.15, true_trend = -0.04)
  expect_equal(mean(d$mean_post_control - d$mean_pre_control), -0.04,
               tolerance = TOL)
})

test_that("baselines, within-study SDs and rho come back as specified", {
  d <- grid_summary(n_did = S_STUDIES, n_rct = 0L, baseline_mean = 0.45,
                    true_effect = -0.15, true_trend = -0.04)
  expect_equal(mean(d$mean_pre_control), 0.45, tolerance = TOL)
  for (col in c("sd_pre_control", "sd_post_control",
                "sd_pre_treatment", "sd_post_treatment")) {
    expect_equal(mean(d[[col]]), 0.12, tolerance = 0.005, info = col)
  }
  expect_equal(mean(d$rho), 0.5, tolerance = 0.02)
  expect_true(all(d$n_control == 200L) && all(d$n_treatment == 200L))
})

test_that("the post-only gap carries theta plus gamma, not theta alone", {
  # A post-only study cannot separate them -- that is the identification
  # problem the whole randomisation model exists to handle -- so the DGP must
  # genuinely confound them.
  d <- grid_summary(n_did = 0L, n_rct = S_STUDIES, true_effect = -0.15,
                    rct_gamma_mean = 0.06, rct_gamma_sd = 0.01)
  expect_equal(mean(d$mean_post_treatment - d$mean_post_control),
               -0.15 + 0.06, tolerance = TOL)
})

test_that("the pre-post change carries theta plus the time trend", {
  d <- grid_summary(n_did = 0L, n_rct = 0L, n_pp = S_STUDIES,
                    true_effect = -0.15, true_trend = -0.04)
  expect_equal(mean(d$mean_post_treatment - d$mean_pre_treatment),
               -0.15 + -0.04, tolerance = TOL)
})

# --- simulate_from_metadid (categories A-W: the bulk of the study) ---------

test_that("the standard simulator recovers theta, beta and the baseline", {
  # metadid::simulate_meta_did() drives every category except X, and nothing
  # here tested it at all before.
  cfg <- scenario("t", dgp = as.list(exact_dgp(n_did = S_STUDIES, n_rct = 0L, n_pp = 0L,
                                               true_effect = -0.15, true_trend = -0.04,
                                               baseline_mean = 0.45)))
  d <- simulate_scenario("t", 303L, cfg)$data$summary_data

  dd <- (d$mean_post_treatment - d$mean_pre_treatment) -
        (d$mean_post_control   - d$mean_pre_control)
  expect_equal(mean(dd), -0.15, tolerance = TOL)
  expect_equal(mean(d$mean_post_control - d$mean_pre_control), -0.04, tolerance = TOL)
  expect_equal(mean(d$mean_pre_control), 0.45, tolerance = TOL)
  expect_equal(mean(d$rho), 0.5, tolerance = 0.02)
})

test_that("the standard simulator generates no baseline imbalance", {
  # Categories A-W assume gamma == 0. If that ever stopped being true their
  # results would shift with no scenario having changed.
  cfg <- scenario("t", dgp = as.list(exact_dgp(n_did = S_STUDIES, n_rct = 0L, n_pp = 0L,
                                               true_effect = -0.15, true_trend = -0.04)))
  d <- simulate_scenario("t", 404L, cfg)$data$summary_data
  expect_equal(mean(d$mean_pre_treatment - d$mean_pre_control), 0, tolerance = TOL)
})

test_that("heterogeneity shows up as between-study spread of the right size", {
  # With sigma_effect switched back on, the spread of per-study double
  # differences must exceed sampling noise by about sigma_effect.
  cfg <- scenario("t", dgp = as.list(mk_dgp(
    n_did = 300L, n_rct = 0L, n_pp = 0L, sigma_effect = 0.05, sigma_trend = 0,
    baseline_sd = 0, within_sd = 0.12, rho = 0.5,
    n_control = 400L, n_treatment = 400L,
    true_effect = -0.15, true_trend = -0.04)))
  d <- simulate_scenario("t", 505L, cfg)$data$summary_data
  dd <- (d$mean_post_treatment - d$mean_pre_treatment) -
        (d$mean_post_control   - d$mean_pre_control)
  # var(observed) = sigma_effect^2 + sampling var; subtract the latter.
  samp_var <- 4 * 0.12^2 * (1 - 0.5) / 400
  expect_equal(sqrt(max(var(dd) - samp_var, 0)), 0.05, tolerance = 0.01)
})

test_that("individual-level data reproduces the summary-level moments", {
  # The two formats must describe the same experiment; a divergence would make
  # the individual-data scenarios silently incomparable with the rest.
  mk <- function(fmt) {
    cfg <- scenario("t",
      dgp = c(list(type = "bespoke", bespoke_fn = "simulate_imbalance_grid"),
              as.list(exact_dgp(n_did = 40L, n_rct = 0L,
                                true_effect = -0.15, true_trend = -0.04))),
      fit = list(data_format = fmt))
    simulate_scenario("t", 606L, cfg)$data
  }
  s <- mk("summary")$summary_data
  i <- mk("individual")$individual_data

  ind_mean <- function(g, tm) {
    i |> filter(group == g, time == tm) |>
      group_by(study_id) |> summarise(m = mean(value), .groups = "drop") |>
      pull(m) |> mean()
  }
  expect_equal(ind_mean("control", "pre"),    mean(s$mean_pre_control),    tolerance = TOL)
  expect_equal(ind_mean("treatment", "post"), mean(s$mean_post_treatment), tolerance = TOL)
})
