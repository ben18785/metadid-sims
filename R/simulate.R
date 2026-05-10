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

# ===========================================================================
# Main dispatcher
# ===========================================================================

simulate_scenario <- function(scenario_id, rep_seed) {

  config <- SCENARIO_CONFIGS[[scenario_id]]
  dgp    <- config$dgp
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

  true_params <- build_true_params(dgp, config$true, fit$normalise_by_baseline)

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
  bm <- dgp$baseline_mean

  params <- tibble(
    treatment_effect_mean_raw        = dgp$true_effect,
    treatment_effect_sd_raw          = dgp$sigma_effect,
    time_trend_mean_raw              = dgp$true_trend,
    time_trend_sd_raw                = dgp$sigma_trend,
    baseline_mean                    = bm,
    treatment_effect_mean_normalised = dgp$true_effect / bm,
    treatment_effect_sd_normalised   = dgp$sigma_effect / bm,
    time_trend_mean_normalised       = dgp$true_trend / bm,
    time_trend_sd_normalised         = dgp$sigma_trend / bm,
    normalised                       = normalised
  )

  # Covariate coefficients
  if (!is.null(dgp$beta_cov) && !is.null(dgp$covariates)) {
    cov_names <- names(dgp$covariates)
    for (i in seq_along(dgp$beta_cov)) {
      params[[paste0("beta_cov_", cov_names[i], "_raw")]] <- dgp$beta_cov[i]
      params[[paste0("beta_cov_", cov_names[i], "_normalised")]] <- dgp$beta_cov[i] / bm
    }
    # True effect at mean covariate value (for centered covariates)
    mean_cov <- colMeans(dgp$covariates)
    effect_at_mean <- dgp$true_effect + sum(mean_cov * dgp$beta_cov)
    params$treatment_effect_at_mean_cov_raw        <- effect_at_mean
    params$treatment_effect_at_mean_cov_normalised  <- effect_at_mean / bm
  }

  # Design-specific effects (F7 etc.)
  if (!is.null(extra_true$delta_rct)) {
    rct_effect <- dgp$true_effect + extra_true$delta_rct
    pp_effect  <- dgp$true_effect + extra_true$delta_pp
    params$delta_rct_raw                          <- extra_true$delta_rct
    params$delta_pp_raw                           <- extra_true$delta_pp
    params$treatment_effect_mean_rct_raw          <- rct_effect
    params$treatment_effect_mean_pp_raw           <- pp_effect
    params$treatment_effect_mean_rct_normalised   <- rct_effect / bm
    params$treatment_effect_mean_pp_normalised    <- pp_effect / bm
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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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

  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)
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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

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

  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)
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

  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)
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

  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)
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
  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)

  list(data = list(summary_data = summary_df), true_params = true_params)
}

# ===========================================================================
# Bespoke DGPs: Edge cases (Category E)
# ===========================================================================

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

  true_params <- build_true_params(dgp, config$true, config$fit$normalise_by_baseline)
  list(data = list(summary_data = summary_df), true_params = true_params)
}
