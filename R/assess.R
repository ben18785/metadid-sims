# Assessment functions for the metadid validation study.
#
# assess_one()          — per-replication: compare posteriors to truth
# aggregate_scenario()  — across replications: coverage, bias, RMSE

library(tibble)
library(dplyr)

# ===========================================================================
# Per-replication assessment
# ===========================================================================

#' Assess one model fit against true parameter values
#'
#' @param posteriors Tibble from extract_posteriors() with columns:
#'   parameter, mean, sd, lo, hi, max_rhat, min_ess_bulk
#' @param true_params One-row tibble from build_true_params()
#' @param fit_config Fit configuration list (to determine which scale to use)
#' @return Tibble with one row per assessed parameter
assess_one <- function(posteriors, true_params, fit_config) {

  normalised <- fit_config$normalise_by_baseline
  suffix <- if (normalised) "_normalised" else "_raw"

  # Map from posterior parameter names to true_params column names
  param_map <- list(
    treatment_effect_mean = paste0("treatment_effect_mean", suffix),
    treatment_effect_sd   = paste0("treatment_effect_sd", suffix),
    time_trend_mean       = paste0("time_trend_mean", suffix),
    time_trend_sd         = paste0("time_trend_sd", suffix)
  )

  # Unnormalised model also estimates baseline

  if (!normalised) {
    param_map$baseline_control_mean <- "baseline_mean"
  }

  # Covariates
  cov_cols <- grep("^beta_cov_.*_raw$", names(true_params), value = TRUE)
  if (length(cov_cols) > 0) {
    cov_names_raw <- sub("_raw$", "", cov_cols)
    for (cn in cov_names_raw) {
      # Posterior parameter name is e.g. "beta_cov[dose]"
      var_name <- sub("beta_cov_", "", cn)
      posterior_name <- paste0("beta_cov[", var_name, "]")
      param_map[[posterior_name]] <- paste0(cn, suffix)
    }

    # If covariates are centered, treatment_effect_mean is effect at mean covariate
    if ("treatment_effect_at_mean_cov_normalised" %in% names(true_params)) {
      param_map$treatment_effect_mean <- paste0("treatment_effect_at_mean_cov", suffix)
    }
  }

  # Design effects
  if ("treatment_effect_mean_rct_normalised" %in% names(true_params) ||
      "treatment_effect_mean_rct_raw" %in% names(true_params)) {
    param_map$treatment_effect_mean_rct <- paste0("treatment_effect_mean_rct", suffix)
    param_map$treatment_effect_mean_pp  <- paste0("treatment_effect_mean_pp", suffix)
  }

  # Assess each mapped parameter
  results <- purrr::imap_dfr(param_map, function(true_col, posterior_name) {
    if (!true_col %in% names(true_params)) return(NULL)
    if (!posterior_name %in% posteriors$parameter) return(NULL)

    true_val <- true_params[[true_col]]
    post_row <- posteriors[posteriors$parameter == posterior_name, ]

    tibble(
      parameter       = posterior_name,
      true_value      = true_val,
      posterior_mean   = post_row$mean,
      posterior_sd     = post_row$sd,
      ci_lo           = post_row$lo,
      ci_hi           = post_row$hi,
      covers          = !is.na(post_row$lo) & post_row$lo < true_val & post_row$hi > true_val,
      bias            = post_row$mean - true_val,
      ci_width        = post_row$hi - post_row$lo,
      max_rhat        = post_row$max_rhat,
      min_ess_bulk    = post_row$min_ess_bulk
    )
  })

  results
}

# ===========================================================================
# Aggregation across replications
# ===========================================================================

#' Aggregate assessment results across replications
#'
#' @param results Tibble of stacked assess_one() outputs across N_REPS,
#'   with columns: scenario_id, model_label, parameter, true_value,
#'   posterior_mean, bias, covers, ci_width, max_rhat, min_ess_bulk
#' @return Summary tibble with one row per (scenario_id, model_label, parameter)
aggregate_scenario <- function(results) {
  results |>
    group_by(scenario_id, model_label, parameter, true_value) |>
    summarise(
      n_reps             = n(),
      empirical_coverage = mean(covers, na.rm = TRUE),
      mean_bias          = mean(bias, na.rm = TRUE),
      median_bias        = median(bias, na.rm = TRUE),
      rmse               = sqrt(mean(bias^2, na.rm = TRUE)),
      mean_ci_width      = mean(ci_width, na.rm = TRUE),
      mean_posterior_sd  = mean(posterior_sd, na.rm = TRUE),
      # Convergence diagnostics
      max_rhat_worst     = max(max_rhat, na.rm = TRUE),
      min_ess_worst      = min(min_ess_bulk, na.rm = TRUE),
      # Fraction with Rhat > 1.05 (problematic fits)
      pct_rhat_bad       = mean(max_rhat > 1.05, na.rm = TRUE),
      .groups = "drop"
    )
}
