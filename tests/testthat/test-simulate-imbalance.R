# simulate_imbalance_grid() must generate what the scenario asked for. If the
# DGP quietly differs from its configuration, the resulting "model failure" is
# really a simulator bug, and no amount of staring at the posterior reveals it.
#
# Large n and a fixed seed keep these deterministic without fitting anything.

sim_summary <- function(..., seed = 101L) {
  cfg <- scenario("t", dgp = c(list(type = "bespoke",
                                    bespoke_fn = "simulate_imbalance_grid"),
                              list(...)))
  simulate_scenario("t", seed, cfg)$data$summary_data
}

test_that("the realised DiD imbalance matches the requested population", {
  d <- sim_summary(n_did = 400L, n_rct = 0L, baseline_sd = 0,
                   did_gamma_mean = 0.08, did_gamma_sd = 0.02)
  realised <- d$mean_pre_treatment - d$mean_pre_control
  expect_equal(mean(realised), 0.08, tolerance = 0.01)
})

test_that("sign flipping preserves magnitude and destroys direction", {
  d <- sim_summary(n_did = 400L, n_rct = 0L, baseline_sd = 0,
                   did_gamma_mean = 0.08, did_gamma_sd = 0.02,
                   did_gamma_sign_flip = TRUE)
  realised <- d$mean_pre_treatment - d$mean_pre_control
  expect_equal(mean(realised), 0, tolerance = 0.015)          # direction gone
  expect_gt(sd(realised), 0.06)                                # magnitude kept
})

test_that("randomised studies get a near-zero imbalance scaled by kappa_true", {
  common <- list(n_did = 300L, n_rct = 0L, n_randomised_did = 300L,
                 baseline_sd = 0, within_sd = 0.12,
                 n_control = 100L, n_treatment = 100L)
  d0 <- do.call(sim_summary, c(common, list(kappa_true = 0)))
  d2 <- do.call(sim_summary, c(common, list(kappa_true = 2)))
  spread <- function(d) sd(d$mean_pre_treatment - d$mean_pre_control)
  # kappa_true = 0 leaves only sampling noise; kappa_true = 2 adds population
  # imbalance of twice the sampling SD on top, so the spread must grow.
  expect_gt(spread(d2), spread(d0) * 1.5)
  expect_equal(mean(d2$mean_pre_treatment - d2$mean_pre_control), 0,
               tolerance = 0.01)
})

test_that("randomisation labels follow n_randomised_*, not the design", {
  d <- sim_summary(n_did = 10L, n_rct = 10L,
                   n_randomised_did = 4L, n_randomised_rct = 10L)
  did <- d[d$design == "did", ]; rct <- d[d$design == "rct", ]
  expect_equal(sum(did$randomisation != "none"), 4L)
  expect_equal(sum(rct$randomisation != "none"), 10L)
  # A randomised DiD and an unrandomised post-only study in one dataset: the
  # combination the design label alone cannot express.
  expect_true(any(did$randomisation == "individual"))
  expect_true(any(did$randomisation == "none"))
})

test_that("mislabel_rate mislabels the stated fraction and nothing else", {
  # Labelled randomised, but drawn from the non-randomised population. The
  # label count must not change -- only the truth behind it.
  d <- sim_summary(n_did = 0L, n_rct = 200L, n_randomised_rct = 200L,
                   baseline_sd = 0, mislabel_rate = 0.5,
                   rct_gamma_mean = 0.10, rct_gamma_sd = 0.01,
                   kappa_true = 0)
  expect_equal(sum(d$randomisation != "none"), 200L)   # all still labelled
  gap <- d$mean_post_treatment - d$mean_post_control
  # Half carry gamma ~ 0.10 on top of the treatment effect, half carry none,
  # so the spread across studies must be far larger than sampling alone.
  expect_gt(sd(gap), 0.03)
})

test_that("the label override changes only the label, not the data", {
  args <- list(n_did = 6L, n_rct = 0L, n_randomised_did = 6L,
               cluster_size = 51, icc = 0.02, kappa_true = 0)
  plain <- do.call(sim_summary, args)
  ovr   <- do.call(sim_summary, c(args, list(randomisation_label_override = "individual")))
  expect_equal(unique(plain$randomisation), "cluster")
  expect_equal(unique(ovr$randomisation),   "individual")
  expect_equal(plain$mean_pre_control, ovr$mean_pre_control)   # identical data
  expect_equal(plain$mean_post_treatment, ovr$mean_post_treatment)
})

test_that("cluster allocation inflates the arm-mean variance by the design effect", {
  # The whole point of the cluster scenarios: arm means must genuinely carry
  # DEFF = 1 + (m-1)*ICC, otherwise the DEFF correction has nothing to correct.
  common <- list(n_did = 250L, n_rct = 0L, n_randomised_did = 250L,
                 baseline_sd = 0, within_sd = 0.12, kappa_true = 0,
                 n_control = 102L, n_treatment = 102L)
  flat <- do.call(sim_summary, common)
  clus <- do.call(sim_summary, c(common, list(cluster_size = 51, icc = 0.02)))
  v <- function(d) var(d$mean_pre_control)
  deff <- 1 + (51 - 1) * 0.02           # = 2
  expect_gt(v(clus) / v(flat), deff * 0.5)   # generous: 250 studies, one draw
  # Marginal individual-level SD is unchanged -- clustering moves variance into
  # the mean, it does not inflate the reported SD.
  expect_equal(mean(clus$sd_pre_control), mean(flat$sd_pre_control),
               tolerance = 0.02)
})

test_that("pre-post studies start at the latent control baseline plus gamma", {
  d <- sim_summary(n_did = 0L, n_rct = 0L, n_pp = 300L, baseline_sd = 0,
                   baseline_mean = 0.45, pp_gamma_mean = 0.08, pp_gamma_sd = 0.01)
  expect_equal(mean(d$mean_pre_treatment), 0.45 + 0.08, tolerance = 0.01)
})

# --- The post-only time trend ------------------------------------------------
# These three guard the defect that produced the "post-only residual bias" in
# X1-X9 and X13-X18: the RCT branch never drew `betas`, so a post-only study's
# arms sat at the latent PRE-treatment level. Its normalised contrast was then
# (gamma + theta)/b, already on the alpha scale that DiD studies use, and
# metadid's alpha + beta -> alpha conversion multiplied in a spurious
# (1 + beta_tilde). The pooled effect came out biased by theta_tilde *
# beta_tilde -- about +0.02 at the default trend, flipping sign with it, and
# invisible to every diagnostic aimed at the model rather than the generator.

test_that("a post-only study's arms are observed after the trend has acted", {
  d <- sim_summary(n_did = 0L, n_rct = 400L, baseline_sd = 0, sigma_trend = 0,
                   baseline_mean = 0.45, true_trend = -0.04)
  d <- d[d$design == "rct", ]
  # The control arm has had no treatment but has had the trend.
  expect_equal(mean(d$mean_post_control), 0.45 - 0.04, tolerance = 0.005)
})

test_that("a post-only study's arm contrast is free of the trend", {
  # gamma == 0 for every study, so the contrast must be theta alone whatever
  # the trend does. This is the complement of the test above: the trend has to
  # be present in the LEVELS and absent from the CONTRAST.
  contrast <- function(trend) {
    d <- sim_summary(n_did = 0L, n_rct = 400L, baseline_sd = 0, sigma_trend = 0,
                     n_randomised_rct = 400L, kappa_true = 0,
                     true_effect = -0.15, sigma_effect = 0, true_trend = trend)
    d <- d[d$design == "rct", ]
    mean(d$mean_post_treatment - d$mean_post_control)
  }
  expect_equal(contrast(-0.04), -0.15, tolerance = 0.005)
  expect_equal(contrast(+0.04), -0.15, tolerance = 0.005)
})

test_that("post-only and DiD studies share the post-period scale", {
  # The invariant whose violation caused the bias: a post-only study's control
  # mean estimates the same quantity as a DiD study's control POST mean. If
  # these drift apart, the two designs are being scored on different scales and
  # any cross-design pooling inherits the difference as bias.
  common <- list(baseline_sd = 0, sigma_trend = 0, baseline_mean = 0.45,
                 true_trend = -0.04)
  dd <- do.call(sim_summary, c(common, list(n_did = 400L, n_rct = 0L)))
  rr <- do.call(sim_summary, c(common, list(n_did = 0L,   n_rct = 400L)))
  expect_equal(mean(rr$mean_post_control[rr$design == "rct"]),
               mean(dd$mean_post_control[dd$design == "did"]),
               tolerance = 0.005)
})

test_that("both data formats describe the same studies", {
  mk <- function(fmt) {
    cfg <- scenario("t",
      dgp = list(type = "bespoke", bespoke_fn = "simulate_imbalance_grid",
                 n_did = 4L, n_rct = 3L, n_pp = 2L, n_randomised_did = 2L),
      fit = list(data_format = fmt))
    simulate_scenario("t", 7L, cfg)$data
  }
  s <- mk("summary")$summary_data
  i <- mk("individual")$individual_data
  expect_equal(sort(unique(s$study_id)), sort(unique(i$study_id)))
  lab <- function(d) d |> distinct(study_id, randomisation) |> arrange(study_id)
  expect_equal(lab(s)$randomisation, lab(i)$randomisation)
  # Individual data must carry the columns metadid's validator requires.
  expect_true(all(c("group", "time", "value", "subject_id") %in% names(i)))
})
