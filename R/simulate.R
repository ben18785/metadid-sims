# Simulation wrappers for the metadid validation study.
#
# Two families of functions:
#   1. simulate_scenario()       — dispatcher that routes to the right simulator
#   2. simulate_from_metadid()   — wraps metadid::simulate_meta_did()
#   3. simulate_bespoke_*()      — hand-crafted DGPs for outlier/violation/edge cases
#
# All functions return a list with:
#   $data        — named list of data frames ready for meta_did()
#                  e.g. list(summary_data = ...) or list(individual_data = ...)
#   $true_params — tibble of true parameter values on BOTH raw and normalised scales

library(tibble)
library(dplyr)
library(metadid)

# Local %||% fallback so simulate.R is independent of rlang attach state.
# (Native %||% is available in base R >= 4.4 but we keep the fallback for
# older R versions and to match the pattern in fit.R / assess.R.)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ===========================================================================
# Main dispatcher
# ===========================================================================

simulate_scenario <- function(scenario_id, rep_seed, config) {

  dgp <- config$dgp
  set.seed(rep_seed)

  if (dgp$type == "metadid") {
    simulate_from_metadid(config)
  } else {
    fn <- match.fun(dgp$bespoke_fn)
    fn(config)
  }
}

# ===========================================================================
# Standard simulation via metadid::simulate_meta_did()
# ===========================================================================

simulate_from_metadid <- function(config) {
  dgp <- config$dgp
  fit <- config$fit

  n_total <- dgp$n_did + dgp$n_rct + dgp$n_pp

  sim <- metadid::simulate_meta_did(
    n_studies        = n_total,
    n_control        = dgp$n_control,
    n_treatment      = dgp$n_treatment,
    true_effect      = dgp$true_effect,
    sigma_effect     = dgp$sigma_effect,
    true_trend       = dgp$true_trend,
    sigma_trend      = dgp$sigma_trend,
    baseline_mean    = dgp$baseline_mean,
    baseline_sd      = dgp$baseline_sd,
    within_sd        = dgp$within_sd,
    rho              = dgp$rho,
    rho_effect_trend = dgp$rho_effect_trend,
    covariates       = dgp$covariates,
    beta_cov         = dgp$beta_cov
  )

  data <- assemble_data(sim, dgp, fit)

  true_params <- build_true_params(dgp, config$true, fit$normalise)

  list(data = data, true_params = true_params)
}

# ===========================================================================
# Data assembly: extract summary/individual data in the right format
# ===========================================================================

assemble_data <- function(sim, dgp, fit) {
  study_ids <- unique(sim$study_id)
  n_total   <- length(study_ids)

  # Partition studies into designs
  idx <- seq_len(n_total)
  did_ids <- study_ids[idx <= dgp$n_did]
  rct_ids <- study_ids[idx > dgp$n_did & idx <= dgp$n_did + dgp$n_rct]
  pp_ids  <- study_ids[idx > dgp$n_did + dgp$n_rct]

  if (fit$data_format == "individual") {
    frames <- list()
    if (length(did_ids) > 0) {
      did_sim <- sim[sim$study_id %in% did_ids, ]
      attr(did_sim, "true_params") <- attr(sim, "true_params")
      frames$did <- metadid::as_individual_did(did_sim)
    }
    if (length(rct_ids) > 0) {
      rct_sim <- sim[sim$study_id %in% rct_ids, ]
      rct_sim$study_id <- paste0("rct_", rct_sim$study_id)
      attr(rct_sim, "true_params") <- attr(sim, "true_params")
      frames$rct <- metadid::as_individual_rct(rct_sim)
    }
    if (length(pp_ids) > 0) {
      pp_sim <- sim[sim$study_id %in% pp_ids, ]
      pp_sim$study_id <- paste0("pp_", pp_sim$study_id)
      attr(pp_sim, "true_params") <- attr(sim, "true_params")
      frames$pp <- metadid::as_individual_pp(pp_sim)
    }
    individual <- bind_rows(frames)
    return(list(individual_data = individual))
  }

  # Summary-level data (default)
  frames <- list()
  if (length(did_ids) > 0) {
    did_sim <- sim[sim$study_id %in% did_ids, ]
    attr(did_sim, "true_params") <- attr(sim, "true_params")
    did_summary <- metadid::as_summary_did(did_sim)
    if (!fit$provide_rho) did_summary$rho <- NA_real_
    frames$did <- did_summary
  }
  if (length(rct_ids) > 0) {
    rct_sim <- sim[sim$study_id %in% rct_ids, ]
    rct_sim$study_id <- paste0("rct_", rct_sim$study_id)
    attr(rct_sim, "true_params") <- attr(sim, "true_params")
    frames$rct <- metadid::as_summary_rct(rct_sim)
  }
  if (length(pp_ids) > 0) {
    pp_sim <- sim[sim$study_id %in% pp_ids, ]
    pp_sim$study_id <- paste0("pp_", pp_sim$study_id)
    attr(pp_sim, "true_params") <- attr(sim, "true_params")
    pp_summary <- metadid::as_summary_pp(pp_sim)
    if (!fit$provide_rho) pp_summary$rho <- NA_real_
    frames$pp <- pp_summary
  }
  summary_df <- bind_rows(frames)
  list(summary_data = summary_df)
}

# ===========================================================================
# True parameter construction
# ===========================================================================

build_true_params <- function(dgp, extra_true, normalised) {
  bm  <- dgp$baseline_mean
  bsd <- if (!is.null(dgp$baseline_sd)) dgp$baseline_sd else 0

  # --- Normalised (percentage-scale) truths ---------------------------------
  # The per-study normalised effect is theta_i / b_i (a unit-free, percentage-
  # scale effect size); the model pools these across studies, so the population
  # estimand is the MEAN of per-study proportional effects, E[theta_i / b_i],
  # with heterogeneity SD[theta_i / b_i].
  #
  # This is NOT true_effect / baseline_mean (= E[theta] / E[baseline]). When
  # baselines vary across studies the two differ by the between-study baseline
  # CV^2 (Jensen's inequality), and averaging raw absolute effects across studies
  # is meaningless when studies are on different scales. Scoring the per-study
  # model against E[theta] / E[baseline] is what produced the spurious "G8 bias".
  # See ben18785/metadid#39.
  #
  # Moments are computed by the delta method (the leading expansion of the
  # per-study ratio theta/b around the means). A full Monte-Carlo expectation
  # over a Gaussian baseline is unusable here because E[theta/b] and Var(theta/b)
  # are dominated by the (unrealistic) near-zero baseline tail of 1/b; the delta
  # approximation matches the regime the model actually sees and is what the
  # historical raw/bm convention was the zeroth-order version of. With
  # baseline_sd == 0 these reduce exactly to the raw/bm values.
  #   E[theta/b]   ~= (E[theta]/bm) * (1 + CV_b^2)
  #   Var(theta/b) ~= ( sigma_theta^2 + E[theta]^2 * CV_b^2 ) / bm^2
  # inv_b_factor = (1 + CV_b^2) rescales a /bm mean into the E[effect/b] mean.
  cv_b2        <- (bsd / bm)^2
  inv_b_factor <- 1 + cv_b2
  te_sd_normalised <- sqrt(dgp$sigma_effect^2 + dgp$true_effect^2 * cv_b2) / bm
  tt_sd_normalised <- sqrt(dgp$sigma_trend^2  + dgp$true_trend^2  * cv_b2) / bm

  params <- tibble(
    treatment_effect_mean_raw        = dgp$true_effect,
    treatment_effect_sd_raw          = dgp$sigma_effect,
    time_trend_mean_raw              = dgp$true_trend,
    time_trend_sd_raw                = dgp$sigma_trend,
    baseline_mean                    = bm,
    treatment_effect_mean_normalised = (dgp$true_effect / bm) * inv_b_factor,
    treatment_effect_sd_normalised   = te_sd_normalised,
    time_trend_mean_normalised       = (dgp$true_trend / bm) * inv_b_factor,
    time_trend_sd_normalised         = tt_sd_normalised,
    normalised                       = normalised
  )

  # Covariate coefficients
  if (!is.null(dgp$beta_cov) && !is.null(dgp$covariates)) {
    cov_names <- names(dgp$covariates)
    for (i in seq_along(dgp$beta_cov)) {
      params[[paste0("beta_cov_", cov_names[i], "_raw")]] <- dgp$beta_cov[i]
      params[[paste0("beta_cov_", cov_names[i], "_normalised")]] <- (dgp$beta_cov[i] / bm) * inv_b_factor
    }
    # True effect at mean covariate value (for centered covariates)
    mean_cov <- colMeans(dgp$covariates)
    effect_at_mean <- dgp$true_effect + sum(mean_cov * dgp$beta_cov)
    params$treatment_effect_at_mean_cov_raw        <- effect_at_mean
    params$treatment_effect_at_mean_cov_normalised  <- (effect_at_mean / bm) * inv_b_factor
  }

  # Design-specific effects (F7 etc.)
  if (!is.null(extra_true$delta_rct)) {
    rct_effect <- dgp$true_effect + extra_true$delta_rct
    pp_effect  <- dgp$true_effect + extra_true$delta_pp
    params$delta_rct_raw                          <- extra_true$delta_rct
    params$delta_pp_raw                           <- extra_true$delta_pp
    params$treatment_effect_mean_rct_raw          <- rct_effect
    params$treatment_effect_mean_pp_raw           <- pp_effect
    params$treatment_effect_mean_rct_normalised   <- (rct_effect / bm) * inv_b_factor
    params$treatment_effect_mean_pp_normalised    <- (pp_effect / bm) * inv_b_factor
  }

  # Rho for effect-trend correlation
  if (dgp$rho_effect_trend != 0) {
    params$rho_effect_trend <- dgp$rho_effect_trend
  }

  params
}

# ===========================================================================
# Bespoke simulation helpers
# ===========================================================================

# Helper: simulate individual bivariate normal data for one study,
# return summary statistics
simulate_one_study <- function(study_id, theta, beta, baseline,
                               n_control, n_treatment, within_sd, rho) {
  Sigma <- within_sd^2 * matrix(c(1, rho, rho, 1), 2, 2)
  mu_ctrl <- c(baseline, baseline + beta)
  mu_trt  <- c(baseline, baseline + beta + theta)

  ctrl <- MASS::mvrnorm(n_control, mu_ctrl, Sigma)
  trt  <- MASS::mvrnorm(n_treatment, mu_trt, Sigma)

  tibble(
    study_id             = study_id,
    n_control            = n_control,
    n_treatment          = n_treatment,
    mean_pre_control     = mean(ctrl[, 1]),
    mean_post_control    = mean(ctrl[, 2]),
    sd_pre_control       = sd(ctrl[, 1]),
    sd_post_control      = sd(ctrl[, 2]),
    mean_pre_treatment   = mean(trt[, 1]),
    mean_post_treatment  = mean(trt[, 2]),
    sd_pre_treatment     = sd(trt[, 1]),
    sd_post_treatment    = sd(trt[, 2]),
    rho                  = (cor(ctrl[, 1], ctrl[, 2]) + cor(trt[, 1], trt[, 2])) / 2
  )
}

# Helper: simulate did summary data from study-level parameter vectors
simulate_did_summary_from_params <- function(study_params, dgp) {
  purrr::pmap_dfr(study_params, function(study_id, theta, beta, baseline, ...) {
    row <- simulate_one_study(
      study_id    = study_id,
      theta       = theta,
      beta        = beta,
      baseline    = baseline,
      n_control   = dgp$n_control,
      n_treatment = dgp$n_treatment,
      within_sd   = dgp$within_sd,
      rho         = dgp$rho
    )
    row$design <- "did"
    row
  })
}

# Helper: simulate one study with multivariate-t within-study errors
simulate_one_study_mvt <- function(study_id, theta, beta, baseline,
                                   n_control, n_treatment, within_sd, rho, df) {
  Sigma <- within_sd^2 * matrix(c(1, rho, rho, 1), 2, 2)
  mu_ctrl <- c(baseline, baseline + beta)
  mu_trt  <- c(baseline, baseline + beta + theta)

  # mvtnorm::rmvt draws from multivariate t with location delta and scale sigma
  ctrl <- mvtnorm::rmvt(n_control, sigma = Sigma * (df - 2) / df, df = df) +
    matrix(mu_ctrl, n_control, 2, byrow = TRUE)
  trt  <- mvtnorm::rmvt(n_treatment, sigma = Sigma * (df - 2) / df, df = df) +
    matrix(mu_trt, n_treatment, 2, byrow = TRUE)

  tibble(
    study_id             = study_id,
    design               = "did",
    n_control            = n_control,
    n_treatment          = n_treatment,
    mean_pre_control     = mean(ctrl[, 1]),
    mean_post_control    = mean(ctrl[, 2]),
    sd_pre_control       = sd(ctrl[, 1]),
    sd_post_control      = sd(ctrl[, 2]),
    mean_pre_treatment   = mean(trt[, 1]),
    mean_post_treatment  = mean(trt[, 2]),
    sd_pre_treatment     = sd(trt[, 1]),
    sd_post_treatment    = sd(trt[, 2]),
    rho                  = (cor(ctrl[, 1], ctrl[, 2]) + cor(trt[, 1], trt[, 2])) / 2
  )
}

# Helper: extract PP summary from a full study simulation
make_pp_summary <- function(row) {
  tibble(
    study_id            = row$study_id,
    design              = "pp",
    n_treatment         = row$n_treatment,
    mean_pre_treatment  = row$mean_pre_treatment,
    sd_pre_treatment    = row$sd_pre_treatment,
    mean_post_treatment = row$mean_post_treatment,
    sd_post_treatment   = row$sd_post_treatment,
    rho                 = row$rho
  )
}

# Helper: extract RCT summary from a full study simulation
make_rct_summary <- function(row) {
  tibble(
    study_id            = row$study_id,
    design              = "rct",
    n_control           = row$n_control,
    n_treatment         = row$n_treatment,
    mean_post_control   = row$mean_post_control,
    sd_post_control     = row$sd_post_control,
    mean_post_treatment = row$mean_post_treatment,
    sd_post_treatment   = row$sd_post_treatment
  )
}

# ===========================================================================
# Bespoke DGPs: Outlier / heavy-tailed (Category C)
# ===========================================================================

simulate_with_outliers <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  # Shift outlier studies
  outlier_idx <- seq_len(dgp$n_outlier)
  shift <- dgp$outlier_shift * dgp$sigma_effect
  if (!is.null(dgp$outlier_direction) && dgp$outlier_direction == "positive") {
    thetas[outlier_idx] <- thetas[outlier_idx] + shift
  } else {
    signs <- sample(c(-1, 1), dgp$n_outlier, replace = TRUE)
    thetas[outlier_idx] <- thetas[outlier_idx] + signs * shift
  }

  study_params <- tibble(
    study_id = paste0("study_", seq_len(n)),
    theta = thetas, beta = betas, baseline = baselines
  )

  summary_df <- simulate_did_summary_from_params(study_params, dgp)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_t_effects <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  # theta_i ~ location-scale t: mean = true_effect, scale = sigma_effect
  thetas    <- dgp$true_effect + dgp$sigma_effect * rt(n, df = dgp$effect_df)
  betas     <- rnorm(n, dgp$true_trend, dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  study_params <- tibble(
    study_id = paste0("study_", seq_len(n)),
    theta = thetas, beta = betas, baseline = baselines
  )

  summary_df <- simulate_did_summary_from_params(study_params, dgp)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_t_likelihood <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  summary_df <- purrr::pmap_dfr(
    tibble(study_id = paste0("study_", seq_len(n)),
           theta = thetas, beta = betas, baseline = baselines),
    function(study_id, theta, beta, baseline, ...) {
      simulate_one_study_mvt(
        study_id, theta, beta, baseline,
        dgp$n_control, dgp$n_treatment, dgp$within_sd, dgp$rho, dgp$within_df
      )
    }
  )

  true_params <- build_true_params(dgp, config$true, config$fit$normalise)
  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_t_trends <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  # beta_i ~ location-scale t
  betas     <- dgp$true_trend + dgp$sigma_trend * rt(n, df = dgp$trend_df)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  study_params <- tibble(
    study_id = paste0("study_", seq_len(n)),
    theta = thetas, beta = betas, baseline = baselines
  )

  summary_df <- simulate_did_summary_from_params(study_params, dgp)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_design_outliers <- function(config) {
  dgp <- config$dgp

  # DiD studies — normal
  did_thetas    <- rnorm(dgp$n_did, dgp$true_effect, dgp$sigma_effect)
  did_betas     <- rnorm(dgp$n_did, dgp$true_trend,  dgp$sigma_trend)
  did_baselines <- rnorm(dgp$n_did, dgp$baseline_mean, dgp$baseline_sd)

  did_params <- tibble(
    study_id = paste0("study_", seq_len(dgp$n_did)),
    theta = did_thetas, beta = did_betas, baseline = did_baselines
  )

  # PP studies — some are outliers
  pp_thetas    <- rnorm(dgp$n_pp, dgp$true_effect, dgp$sigma_effect)
  pp_betas     <- rnorm(dgp$n_pp, dgp$true_trend,  dgp$sigma_trend)
  pp_baselines <- rnorm(dgp$n_pp, dgp$baseline_mean, dgp$baseline_sd)

  outlier_idx <- seq_len(dgp$n_outlier_pp)
  shift <- dgp$outlier_shift * dgp$sigma_effect
  pp_thetas[outlier_idx] <- pp_thetas[outlier_idx] +
    sample(c(-1, 1), dgp$n_outlier_pp, replace = TRUE) * shift

  pp_params <- tibble(
    study_id = paste0("pp_study_", seq_len(dgp$n_pp)),
    theta = pp_thetas, beta = pp_betas, baseline = pp_baselines
  )

  # Generate summaries
  did_summary <- simulate_did_summary_from_params(did_params, dgp)

  pp_summary <- purrr::pmap_dfr(pp_params, function(study_id, theta, beta, baseline, ...) {
    row <- simulate_one_study(study_id, theta, beta, baseline,
                              dgp$n_control, dgp$n_treatment, dgp$within_sd, dgp$rho)
    make_pp_summary(row)
  })

  summary_df <- bind_rows(did_summary, pp_summary)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

# ===========================================================================
# Bespoke DGPs: Assumption violations (Category D)
# ===========================================================================

simulate_design_offsets <- function(config) {
  dgp <- config$dgp

  # DiD studies: standard effect
  did_thetas    <- rnorm(dgp$n_did, dgp$true_effect, dgp$sigma_effect)
  did_betas     <- rnorm(dgp$n_did, dgp$true_trend,  dgp$sigma_trend)
  did_baselines <- rnorm(dgp$n_did, dgp$baseline_mean, dgp$baseline_sd)

  did_params <- tibble(
    study_id = paste0("study_", seq_len(dgp$n_did)),
    theta = did_thetas, beta = did_betas, baseline = did_baselines
  )
  did_summary <- simulate_did_summary_from_params(did_params, dgp)

  # RCT studies: shifted effect
  rct_true_effect <- dgp$true_effect + dgp$delta_rct
  rct_thetas    <- rnorm(dgp$n_rct, rct_true_effect, dgp$sigma_effect)
  rct_betas     <- rnorm(dgp$n_rct, dgp$true_trend,  dgp$sigma_trend)
  rct_baselines <- rnorm(dgp$n_rct, dgp$baseline_mean, dgp$baseline_sd)

  rct_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("rct_", seq_len(dgp$n_rct)),
           theta = rct_thetas, beta = rct_betas, baseline = rct_baselines),
    function(study_id, theta, beta, baseline, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, dgp$within_sd, dgp$rho)
      make_rct_summary(row)
    }
  )

  # PP studies: shifted effect
  pp_true_effect <- dgp$true_effect + dgp$delta_pp
  pp_thetas    <- rnorm(dgp$n_pp, pp_true_effect, dgp$sigma_effect)
  pp_betas     <- rnorm(dgp$n_pp, dgp$true_trend, dgp$sigma_trend)
  pp_baselines <- rnorm(dgp$n_pp, dgp$baseline_mean, dgp$baseline_sd)

  pp_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("pp_", seq_len(dgp$n_pp)),
           theta = pp_thetas, beta = pp_betas, baseline = pp_baselines),
    function(study_id, theta, beta, baseline, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, dgp$within_sd, dgp$rho)
      make_pp_summary(row)
    }
  )

  summary_df <- bind_rows(did_summary, rct_summary, pp_summary)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_divergent_trends <- function(config) {
  dgp <- config$dgp

  # DiD studies: use did_trend
  did_thetas    <- rnorm(dgp$n_did, dgp$true_effect, dgp$sigma_effect)
  did_betas     <- rnorm(dgp$n_did, dgp$did_trend, dgp$sigma_trend)
  did_baselines <- rnorm(dgp$n_did, dgp$baseline_mean, dgp$baseline_sd)

  did_params <- tibble(
    study_id = paste0("study_", seq_len(dgp$n_did)),
    theta = did_thetas, beta = did_betas, baseline = did_baselines
  )
  did_summary <- simulate_did_summary_from_params(did_params, dgp)

  # PP studies: use pp_trend
  pp_thetas    <- rnorm(dgp$n_pp, dgp$true_effect, dgp$sigma_effect)
  pp_betas     <- rnorm(dgp$n_pp, dgp$pp_trend, dgp$sigma_trend)
  pp_baselines <- rnorm(dgp$n_pp, dgp$baseline_mean, dgp$baseline_sd)

  pp_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("pp_", seq_len(dgp$n_pp)),
           theta = pp_thetas, beta = pp_betas, baseline = pp_baselines),
    function(study_id, theta, beta, baseline, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, dgp$within_sd, dgp$rho)
      make_pp_summary(row)
    }
  )

  summary_df <- bind_rows(did_summary, pp_summary)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_hetero_variance <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)
  within_sds <- rgamma(n, shape = dgp$within_sd_shape, rate = dgp$within_sd_rate)

  summary_df <- purrr::pmap_dfr(
    tibble(study_id = paste0("study_", seq_len(n)),
           theta = thetas, beta = betas, baseline = baselines,
           w_sd = within_sds),
    function(study_id, theta, beta, baseline, w_sd, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, w_sd, dgp$rho)
      row$design <- "did"
      row
    }
  )

  true_params <- build_true_params(dgp, config$true, config$fit$normalise)
  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_varying_rho <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)
  rhos      <- rbeta(n, dgp$rho_alpha, dgp$rho_beta)

  summary_df <- purrr::pmap_dfr(
    tibble(study_id = paste0("study_", seq_len(n)),
           theta = thetas, beta = betas, baseline = baselines,
           rho_i = rhos),
    function(study_id, theta, beta, baseline, rho_i, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, dgp$within_sd, rho_i)
      row$design <- "did"
      row
    }
  )

  true_params <- build_true_params(dgp, config$true, config$fit$normalise)
  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_size_effect_corr <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did

  n_controls   <- sample(30:200, n, replace = TRUE)
  n_treatments <- sample(30:200, n, replace = TRUE)

  mean_n <- mean(c(n_controls, n_treatments)) / 2
  thetas <- rnorm(n, dgp$true_effect + dgp$size_effect_slope * (n_treatments - mean_n),
                  dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend, dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  summary_df <- purrr::pmap_dfr(
    tibble(study_id = paste0("study_", seq_len(n)),
           theta = thetas, beta = betas, baseline = baselines,
           nc = n_controls, nt = n_treatments),
    function(study_id, theta, beta, baseline, nc, nt, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                nc, nt, dgp$within_sd, dgp$rho)
      row$design <- "did"
      row
    }
  )

  true_params <- build_true_params(dgp, config$true, config$fit$normalise)
  list(data = list(summary_data = summary_df), true_params = true_params)
}

simulate_rct_imbalance <- function(config) {
  dgp <- config$dgp

  # DiD studies — standard
  did_thetas    <- rnorm(dgp$n_did, dgp$true_effect, dgp$sigma_effect)
  did_betas     <- rnorm(dgp$n_did, dgp$true_trend,  dgp$sigma_trend)
  did_baselines <- rnorm(dgp$n_did, dgp$baseline_mean, dgp$baseline_sd)

  did_params <- tibble(
    study_id = paste0("study_", seq_len(dgp$n_did)),
    theta = did_thetas, beta = did_betas, baseline = did_baselines
  )
  did_summary <- simulate_did_summary_from_params(did_params, dgp)

  # RCT studies — with baseline imbalance (gamma != 0)
  rct_thetas    <- rnorm(dgp$n_rct, dgp$true_effect, dgp$sigma_effect)
  rct_baselines <- rnorm(dgp$n_rct, dgp$baseline_mean, dgp$baseline_sd)
  rct_gammas    <- rnorm(dgp$n_rct, dgp$rct_gamma_mean, dgp$rct_gamma_sd)

  rct_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("rct_", seq_len(dgp$n_rct)),
           theta = rct_thetas, baseline = rct_baselines, gamma = rct_gammas),
    function(study_id, theta, baseline, gamma, ...) {
      # Cross-sectional: control post ~ N(baseline, sigma^2),
      #                  treatment post ~ N(baseline + gamma + theta, sigma^2)
      ctrl_post <- rnorm(dgp$n_control, baseline, dgp$within_sd)
      trt_post  <- rnorm(dgp$n_treatment, baseline + gamma + theta, dgp$within_sd)
      tibble(
        study_id            = study_id,
        design              = "rct",
        n_control           = dgp$n_control,
        n_treatment         = dgp$n_treatment,
        mean_post_control   = mean(ctrl_post),
        sd_post_control     = sd(ctrl_post),
        mean_post_treatment = mean(trt_post),
        sd_post_treatment   = sd(trt_post)
      )
    }
  )

  summary_df <- bind_rows(did_summary, rct_summary)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

# Simulate DiD + RCT where DiD studies have larger within-study baseline
# imbalance (treatment arm starts at baseline + gamma) than RCT studies.
# DiD differencing recovers theta regardless of gamma, but a naive model that
# pools imbalance across designs may incorrectly inflate RCT imbalance
# estimates when DiD imbalance is large.
#
# DGP parameters (beyond defaults):
#   did_gamma_mean / did_gamma_sd  — imbalance distribution for DiD
#   rct_gamma_mean / rct_gamma_sd  — imbalance distribution for RCT
simulate_did_rct_imbalance <- function(config) {
  dgp <- config$dgp

  # DiD studies: treatment arm starts at baseline + gamma
  did_thetas    <- rnorm(dgp$n_did, dgp$true_effect, dgp$sigma_effect)
  did_betas     <- rnorm(dgp$n_did, dgp$true_trend,  dgp$sigma_trend)
  did_baselines <- rnorm(dgp$n_did, dgp$baseline_mean, dgp$baseline_sd)
  did_gammas    <- rnorm(dgp$n_did, dgp$did_gamma_mean, dgp$did_gamma_sd)

  did_summary <- purrr::pmap_dfr(
    tibble(
      study_id = paste0("did_", seq_len(dgp$n_did)),
      theta = did_thetas, beta = did_betas,
      baseline = did_baselines, gamma = did_gammas
    ),
    function(study_id, theta, beta, baseline, gamma, ...) {
      Sigma   <- dgp$within_sd^2 * matrix(c(1, dgp$rho, dgp$rho, 1), 2, 2)
      mu_ctrl <- c(baseline,         baseline + beta)
      mu_trt  <- c(baseline + gamma, baseline + gamma + beta + theta)
      ctrl <- MASS::mvrnorm(dgp$n_control,   mu_ctrl, Sigma)
      trt  <- MASS::mvrnorm(dgp$n_treatment, mu_trt,  Sigma)
      tibble(
        study_id             = study_id,
        design               = "did",
        n_control            = dgp$n_control,
        n_treatment          = dgp$n_treatment,
        mean_pre_control     = mean(ctrl[, 1]),
        mean_post_control    = mean(ctrl[, 2]),
        sd_pre_control       = sd(ctrl[, 1]),
        sd_post_control      = sd(ctrl[, 2]),
        mean_pre_treatment   = mean(trt[, 1]),
        mean_post_treatment  = mean(trt[, 2]),
        sd_pre_treatment     = sd(trt[, 1]),
        sd_post_treatment    = sd(trt[, 2]),
        rho = (cor(ctrl[, 1], ctrl[, 2]) + cor(trt[, 1], trt[, 2])) / 2
      )
    }
  )

  # RCT studies: small/zero imbalance
  rct_thetas    <- rnorm(dgp$n_rct, dgp$true_effect, dgp$sigma_effect)
  rct_baselines <- rnorm(dgp$n_rct, dgp$baseline_mean, dgp$baseline_sd)
  rct_gammas    <- rnorm(dgp$n_rct, dgp$rct_gamma_mean, dgp$rct_gamma_sd)

  rct_summary <- purrr::pmap_dfr(
    tibble(
      study_id = paste0("rct_", seq_len(dgp$n_rct)),
      theta = rct_thetas, baseline = rct_baselines, gamma = rct_gammas
    ),
    function(study_id, theta, baseline, gamma, ...) {
      ctrl_post <- rnorm(dgp$n_control,   baseline,               dgp$within_sd)
      trt_post  <- rnorm(dgp$n_treatment, baseline + gamma + theta, dgp$within_sd)
      tibble(
        study_id            = study_id,
        design              = "rct",
        n_control           = dgp$n_control,
        n_treatment         = dgp$n_treatment,
        mean_post_control   = mean(ctrl_post),
        sd_post_control     = sd(ctrl_post),
        mean_post_treatment = mean(trt_post),
        sd_post_treatment   = sd(trt_post)
      )
    }
  )

  summary_df  <- bind_rows(did_summary, rct_summary)
  true_params <- build_true_params(dgp, config$true, config$fit$normalise)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

# ===========================================================================
# Bespoke DGPs: Edge cases (Category E)
# ===========================================================================

# ===========================================================================
# Bespoke DGPs: Time trend distribution misspecification (Category H)
# ===========================================================================

#' Skewed (log-normal) time trends: beta_i ~ -LogNormal
#'
#' The model assumes beta_i ~ N(mu, sigma^2), but here beta_i is drawn from a
#' negated log-normal with the same mean as true_trend. lognormal_sigma
#' (on the log scale) controls the degree of skewness.
simulate_lognormal_trends <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did + dgp$n_pp

  # Derive log-scale location so that E[beta_i] = true_trend
  lognormal_mu <- log(abs(dgp$true_trend)) - dgp$lognormal_sigma^2 / 2

  thetas    <- rnorm(n, dgp$true_effect,  dgp$sigma_effect)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)
  betas     <- -rlnorm(n, lognormal_mu, dgp$lognormal_sigma)

  did_idx <- seq_len(dgp$n_did)
  pp_idx  <- seq(dgp$n_did + 1L, n)

  did_summary <- simulate_did_summary_from_params(
    tibble(study_id = paste0("did_", seq_len(dgp$n_did)),
           theta = thetas[did_idx], beta = betas[did_idx],
           baseline = baselines[did_idx]),
    dgp
  )

  pp_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("pp_", seq_len(dgp$n_pp)),
           theta = thetas[pp_idx], beta = betas[pp_idx],
           baseline = baselines[pp_idx]),
    function(study_id, theta, beta, baseline, ...) {
      make_pp_summary(simulate_one_study(study_id, theta, beta, baseline,
                                         dgp$n_control, dgp$n_treatment,
                                         dgp$within_sd, dgp$rho))
    }
  )

  list(
    data        = list(summary_data = bind_rows(did_summary, pp_summary)),
    true_params = build_true_params(dgp, config$true, config$fit$normalise)
  )
}

#' Bimodal time trends: mixture of two normals
#'
#' beta_i is drawn from a 50/50 mixture of N(trend_low, sigma_within^2) and
#' N(trend_high, sigma_within^2). Both components have the same within-group SD.
#' true_trend in dgp should equal mix_prob*trend_high + (1-mix_prob)*trend_low.
simulate_bimodal_trends <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did + dgp$n_pp

  component <- rbinom(n, 1L, dgp$mix_prob)
  betas     <- ifelse(
    component == 1L,
    rnorm(n, dgp$trend_high, dgp$sigma_within),
    rnorm(n, dgp$trend_low,  dgp$sigma_within)
  )
  thetas    <- rnorm(n, dgp$true_effect,  dgp$sigma_effect)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)

  did_idx <- seq_len(dgp$n_did)
  pp_idx  <- seq(dgp$n_did + 1L, n)

  did_summary <- simulate_did_summary_from_params(
    tibble(study_id = paste0("did_", seq_len(dgp$n_did)),
           theta = thetas[did_idx], beta = betas[did_idx],
           baseline = baselines[did_idx]),
    dgp
  )

  pp_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("pp_", seq_len(dgp$n_pp)),
           theta = thetas[pp_idx], beta = betas[pp_idx],
           baseline = baselines[pp_idx]),
    function(study_id, theta, beta, baseline, ...) {
      make_pp_summary(simulate_one_study(study_id, theta, beta, baseline,
                                         dgp$n_control, dgp$n_treatment,
                                         dgp$within_sd, dgp$rho))
    }
  )

  list(
    data        = list(summary_data = bind_rows(did_summary, pp_summary)),
    true_params = build_true_params(dgp, config$true, config$fit$normalise)
  )
}

#' Trend-size confounding: DiD studies have small n, PP studies have large n
#'
#' Trend is negatively correlated with study size (larger n -> trend closer to
#' zero). DiD studies are drawn with n_treatment ~ U(n_small_min, n_small_max)
#' and PP studies with n_treatment ~ U(n_large_min, n_large_max). This creates
#' a design-confounded trend: the model learns trend primarily from small DiD
#' studies and misapplies it to the larger PP studies.
simulate_trend_size_corr <- function(config) {
  dgp <- config$dgp

  n_did <- dgp$n_did
  n_pp  <- dgp$n_pp

  # Study-specific sample sizes
  n_did_vec <- sample(dgp$n_small_min:dgp$n_small_max, n_did, replace = TRUE)
  n_pp_vec  <- sample(dgp$n_large_min:dgp$n_large_max, n_pp,  replace = TRUE)
  all_n     <- c(n_did_vec, n_pp_vec)
  mean_n    <- mean(all_n)

  # Trend proportional to deviation from mean study size
  betas_did <- dgp$true_trend + dgp$trend_size_slope * (n_did_vec - mean_n) +
               rnorm(n_did, 0, dgp$sigma_trend)
  betas_pp  <- dgp$true_trend + dgp$trend_size_slope * (n_pp_vec  - mean_n) +
               rnorm(n_pp,  0, dgp$sigma_trend)

  thetas_did    <- rnorm(n_did, dgp$true_effect,  dgp$sigma_effect)
  thetas_pp     <- rnorm(n_pp,  dgp$true_effect,  dgp$sigma_effect)
  baselines_did <- rnorm(n_did, dgp$baseline_mean, dgp$baseline_sd)
  baselines_pp  <- rnorm(n_pp,  dgp$baseline_mean, dgp$baseline_sd)

  did_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("did_", seq_len(n_did)),
           theta = thetas_did, beta = betas_did,
           baseline = baselines_did, nc = n_did_vec),
    function(study_id, theta, beta, baseline, nc, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                nc, nc, dgp$within_sd, dgp$rho)
      row$design <- "did"
      row
    }
  )

  pp_summary <- purrr::pmap_dfr(
    tibble(study_id = paste0("pp_", seq_len(n_pp)),
           theta = thetas_pp, beta = betas_pp,
           baseline = baselines_pp, nc = n_pp_vec),
    function(study_id, theta, beta, baseline, nc, ...) {
      make_pp_summary(simulate_one_study(study_id, theta, beta, baseline,
                                         nc, nc, dgp$within_sd, dgp$rho))
    }
  )

  list(
    data        = list(summary_data = bind_rows(did_summary, pp_summary)),
    true_params = build_true_params(dgp, config$true, config$fit$normalise)
  )
}

simulate_extreme_rho <- function(config) {
  dgp <- config$dgp
  n   <- dgp$n_did
  n_half <- n %/% 2

  thetas    <- rnorm(n, dgp$true_effect, dgp$sigma_effect)
  betas     <- rnorm(n, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n, dgp$baseline_mean, dgp$baseline_sd)
  rhos      <- c(rep(dgp$rho_values[1], n_half),
                 rep(dgp$rho_values[2], n - n_half))

  summary_df <- purrr::pmap_dfr(
    tibble(study_id = paste0("study_", seq_len(n)),
           theta = thetas, beta = betas, baseline = baselines,
           rho_i = rhos),
    function(study_id, theta, beta, baseline, rho_i, ...) {
      row <- simulate_one_study(study_id, theta, beta, baseline,
                                dgp$n_control, dgp$n_treatment, dgp$within_sd, rho_i)
      row$design <- "did"
      row
    }
  )

  true_params <- build_true_params(dgp, config$true, config$fit$normalise)
  list(data = list(summary_data = summary_df), true_params = true_params)
}

# ===========================================================================
# Bespoke DGPs: Multiplicative covariate (Category I)
# ===========================================================================

# Generate raw bivariate observations (pre, post) for one study × one arm.
.gen_arm_observations <- function(n, mu_pre, mu_post, within_sd, rho) {
  Sigma <- within_sd^2 * matrix(c(1, rho, rho, 1), 2, 2)
  MASS::mvrnorm(n, c(mu_pre, mu_post), Sigma)
}

# Long-form rows for one study, given pre-computed observation matrices.
# Used by the individual-level path of simulate_multiplicative_levels.
.long_rows_one_study <- function(study_id, design, ctrl_obs, trt_obs) {
  n_c <- nrow(ctrl_obs %||% matrix(0, 0, 2))
  n_t <- nrow(trt_obs)

  ctrl_long <- if (n_c > 0) {
    tibble(
      study_id   = study_id,
      subject_id = paste0(study_id, "_c", rep(seq_len(n_c), times = 2)),
      design     = design,
      group      = "control",
      time       = rep(c("pre", "post"), each = n_c),
      value      = c(ctrl_obs[, 1], ctrl_obs[, 2])
    )
  } else tibble()

  trt_long <- tibble(
    study_id   = study_id,
    subject_id = paste0(study_id, "_t", rep(seq_len(n_t), times = 2)),
    design     = design,
    group      = "treatment",
    time       = rep(c("pre", "post"), each = n_t),
    value      = c(trt_obs[, 1], trt_obs[, 2])
  )

  bind_rows(ctrl_long, trt_long)
}

# Simulator for the I-category multiplicative-covariate scenarios.
#
# Per-study true effect θ_i is drawn from
#   N(m_{level(i)} · (true_effect + X_cov_i · beta_cov), sigma_effect²)
# where m_k is the per-level multiplier (m_1 = 1 by convention; reference
# level).
#
# Required DGP fields:
#   level_assignments — integer vector of length n_did + n_rct + n_pp,
#                       with values in 1:length(level_multipliers).
#                       Ordering matches design ordering (DiD, then RCT,
#                       then PP), matching the convention used elsewhere
#                       in this file.
#   level_multipliers — numeric vector of length K. The first entry is
#                       the reference level and MUST equal 1.
# Optional DGP fields:
#   level_names       — character vector of length K, used for the names
#                       of the effect_multiplier_<name> truth columns and
#                       the level column emitted in the data. Defaults to
#                       letters[1:K] so e.g. "a" is the reference.
#   level2_assignments / level2_multipliers / level2_names
#                     — same as the level_* fields but for an optional
#                       SECOND multiplicative covariate (emitted as a
#                       `level2` column). When present, each study's
#                       overall factor is the PRODUCT of its two
#                       per-covariate multipliers, matching metadid's
#                       `multiplicative_covariate = ~ level + level2`
#                       product structure. With two covariates the truth
#                       columns are named effect_multiplier_<column>:<name>
#                       to match metadid's two-covariate posterior labels
#                       effect_multiplier[<column>:<level>].
#   covariates / beta_cov — standard additive covariate fields (passed
#                       through as in simulate_meta_did).
#
# Output shape mirrors the other bespoke simulators: list(data, true_params).
# data contains either summary_data or individual_data depending on
# config$fit$data_format. Both carry a `level` column that is constant within
# each study; if covariates are present, the covariate columns are carried
# through too.
#
# Used by: I1 (binary, balanced), I2 (three-level), I3 (multiplicative +
# additive), I4 (mixed designs), I5 (omit-multiplier comparator), I6
# (spurious multiplier, all studies at reference), I7 (individual-level
# twin of I1), I8 (two-covariate product, binary × binary), I9
# (two-covariate product, binary × three-level).
simulate_multiplicative_levels <- function(config) {
  dgp <- config$dgp
  fit <- config$fit

  n_did <- dgp$n_did %||% 0L
  n_rct <- dgp$n_rct %||% 0L
  n_pp  <- dgp$n_pp  %||% 0L
  n_total <- n_did + n_rct + n_pp

  level_assignments <- as.integer(dgp$level_assignments)
  multipliers       <- dgp$level_multipliers
  K                 <- length(multipliers)
  level_names       <- dgp$level_names %||% letters[seq_len(K)]

  stopifnot(
    length(level_assignments) == n_total,
    all(level_assignments >= 1L & level_assignments <= K),
    abs(multipliers[1] - 1) < 1e-9,        # reference must be 1
    length(level_names) == K
  )

  # Optional second multiplicative covariate (product structure)
  has_second <- !is.null(dgp$level2_multipliers)
  if (has_second) {
    level2_assignments <- as.integer(dgp$level2_assignments)
    multipliers2       <- dgp$level2_multipliers
    K2                 <- length(multipliers2)
    level2_names       <- dgp$level2_names %||% letters[seq_len(K2)]
    stopifnot(
      length(level2_assignments) == n_total,
      all(level2_assignments >= 1L & level2_assignments <= K2),
      abs(multipliers2[1] - 1) < 1e-9,     # reference must be 1
      length(level2_names) == K2
    )
    per_study_level2_name <- level2_names[level2_assignments]
  }

  per_study_multiplier <- multipliers[level_assignments]
  if (has_second) {
    per_study_multiplier <- per_study_multiplier * multipliers2[level2_assignments]
  }
  per_study_level_name <- level_names[level_assignments]

  # Optional additive covariate contribution to per-study true effect
  if (!is.null(dgp$covariates) && !is.null(dgp$beta_cov)) {
    X <- as.matrix(dgp$covariates)
    cov_effect <- as.numeric(X %*% as.numeric(dgp$beta_cov))
  } else {
    cov_effect <- rep(0, n_total)
  }

  # Per-study true effect mean: m_i · (μ_θ + X_cov_i · β_cov)
  multiplied_mean <- per_study_multiplier * (dgp$true_effect + cov_effect)

  thetas    <- rnorm(n_total, multiplied_mean, dgp$sigma_effect)
  betas     <- rnorm(n_total, dgp$true_trend,  dgp$sigma_trend)
  baselines <- rnorm(n_total, dgp$baseline_mean, dgp$baseline_sd)

  # Study IDs and design assignments (DiD first, then RCT, then PP)
  ids_did <- if (n_did > 0) paste0("did_", seq_len(n_did)) else character()
  ids_rct <- if (n_rct > 0) paste0("rct_", seq_len(n_rct)) else character()
  ids_pp  <- if (n_pp  > 0) paste0("pp_",  seq_len(n_pp))  else character()
  ids     <- c(ids_did, ids_rct, ids_pp)
  designs <- c(rep("did", n_did), rep("rct", n_rct), rep("pp", n_pp))

  # Generate per-study observation matrices once; reuse for both summary
  # and individual output paths.
  obs_list <- vector("list", n_total)
  for (i in seq_len(n_total)) {
    mu_ctrl_pre  <- baselines[i]
    mu_ctrl_post <- baselines[i] + betas[i]
    mu_trt_pre   <- baselines[i]
    mu_trt_post  <- baselines[i] + betas[i] + thetas[i]

    ctrl_obs <- .gen_arm_observations(dgp$n_control,
                                       mu_ctrl_pre, mu_ctrl_post,
                                       dgp$within_sd, dgp$rho)
    trt_obs  <- .gen_arm_observations(dgp$n_treatment,
                                       mu_trt_pre,  mu_trt_post,
                                       dgp$within_sd, dgp$rho)
    obs_list[[i]] <- list(ctrl = ctrl_obs, trt = trt_obs)
  }

  if (fit$data_format == "individual") {
    long_frames <- list()
    for (i in seq_len(n_total)) {
      design <- designs[i]
      ctrl_obs <- obs_list[[i]]$ctrl
      trt_obs  <- obs_list[[i]]$trt

      if (design == "rct") {
        # RCT: post-period observations only, both arms
        n_c <- nrow(ctrl_obs); n_t <- nrow(trt_obs)
        rows <- tibble(
          study_id   = ids[i],
          design     = "rct",
          group      = c(rep("control", n_c), rep("treatment", n_t)),
          time       = "post",
          value      = c(ctrl_obs[, 2], trt_obs[, 2])
        )
      } else if (design == "pp") {
        # PP: treatment arm only, both periods
        n_t <- nrow(trt_obs)
        rows <- tibble(
          study_id   = ids[i],
          subject_id = paste0(ids[i], "_t", rep(seq_len(n_t), times = 2)),
          design     = "pp",
          group      = "treatment",
          time       = rep(c("pre", "post"), each = n_t),
          value      = c(trt_obs[, 1], trt_obs[, 2])
        )
      } else {
        # DiD: both arms, both periods
        rows <- .long_rows_one_study(ids[i], "did", ctrl_obs, trt_obs)
      }

      # Attach study-level columns (level and any additive covariates)
      rows$level <- per_study_level_name[i]
      if (has_second) rows$level2 <- per_study_level2_name[i]
      if (!is.null(dgp$covariates)) {
        for (cn in names(dgp$covariates)) {
          rows[[cn]] <- dgp$covariates[[cn]][i]
        }
      }
      long_frames[[i]] <- rows
    }
    individual <- bind_rows(long_frames)
    data <- list(individual_data = individual)

  } else {
    # Summary path
    summary_frames <- list()
    for (i in seq_len(n_total)) {
      design <- designs[i]
      ctrl_obs <- obs_list[[i]]$ctrl
      trt_obs  <- obs_list[[i]]$trt

      if (design == "did") {
        row <- tibble(
          study_id             = ids[i],
          design               = "did",
          n_control            = dgp$n_control,
          n_treatment          = dgp$n_treatment,
          mean_pre_control     = mean(ctrl_obs[, 1]),
          mean_post_control    = mean(ctrl_obs[, 2]),
          sd_pre_control       = sd(ctrl_obs[, 1]),
          sd_post_control      = sd(ctrl_obs[, 2]),
          mean_pre_treatment   = mean(trt_obs[, 1]),
          mean_post_treatment  = mean(trt_obs[, 2]),
          sd_pre_treatment     = sd(trt_obs[, 1]),
          sd_post_treatment    = sd(trt_obs[, 2]),
          rho                  = (cor(ctrl_obs[, 1], ctrl_obs[, 2]) +
                                    cor(trt_obs[, 1], trt_obs[, 2])) / 2
        )
      } else if (design == "rct") {
        row <- tibble(
          study_id            = ids[i],
          design              = "rct",
          n_control           = dgp$n_control,
          n_treatment         = dgp$n_treatment,
          mean_post_control   = mean(ctrl_obs[, 2]),
          sd_post_control     = sd(ctrl_obs[, 2]),
          mean_post_treatment = mean(trt_obs[, 2]),
          sd_post_treatment   = sd(trt_obs[, 2])
        )
      } else {  # pp
        row <- tibble(
          study_id            = ids[i],
          design              = "pp",
          n_treatment         = dgp$n_treatment,
          mean_pre_treatment  = mean(trt_obs[, 1]),
          sd_pre_treatment    = sd(trt_obs[, 1]),
          mean_post_treatment = mean(trt_obs[, 2]),
          sd_post_treatment   = sd(trt_obs[, 2]),
          rho                 = cor(trt_obs[, 1], trt_obs[, 2])
        )
      }
      row$level <- per_study_level_name[i]
      if (has_second) row$level2 <- per_study_level2_name[i]
      if (!is.null(dgp$covariates)) {
        for (cn in names(dgp$covariates)) {
          row[[cn]] <- dgp$covariates[[cn]][i]
        }
      }
      summary_frames[[i]] <- row
    }
    summary_df <- bind_rows(summary_frames)
    if (!fit$provide_rho && "rho" %in% names(summary_df)) {
      summary_df$rho <- NA_real_
    }
    data <- list(summary_data = summary_df)
  }

  # Build truth params, appending one column per level for effect_multiplier.
  # With one covariate metadid labels posteriors effect_multiplier[<level>];
  # with two it labels them effect_multiplier[<column>:<level>], so the truth
  # columns follow the same convention for assess_one()'s mapping.
  true_params <- build_true_params(dgp, config$true, fit$normalise)
  if (has_second) {
    for (k in seq_len(K)) {
      true_params[[paste0("effect_multiplier_level:", level_names[k])]] <- multipliers[k]
    }
    for (k in seq_len(K2)) {
      true_params[[paste0("effect_multiplier_level2:", level2_names[k])]] <- multipliers2[k]
    }
  } else {
    for (k in seq_len(K)) {
      col <- paste0("effect_multiplier_", level_names[k])
      true_params[[col]] <- multipliers[k]
    }
  }

  list(data = data, true_params = true_params)
}

# Helper for scenario configs: build a balanced level_assignments vector
# for n_total studies and K levels (interleaved as evenly as possible).
balanced_level_assignments <- function(n_total, K) {
  rep_len(seq_len(K), n_total)
}
