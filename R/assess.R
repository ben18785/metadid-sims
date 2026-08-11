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

  normalised <- fit_config$normalise
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

  # Multiplicative-covariate factor multipliers (one per level, including
  # the reference). The reference level is fixed at 1 in the fit and isn't
  # sampled, so there's no `effect_multiplier[<reference>]` row in the
  # posterior — the lookup below returns NULL for that mapping and it gets
  # silently dropped. Non-reference levels appear as
  # `effect_multiplier[<level_name>]` in summary().
  #
  # The multiplier is a dimensionless ratio (same value on raw and
  # normalised scales), so no _raw / _normalised suffix is involved.
  mult_cols <- grep("^effect_multiplier_", names(true_params), value = TRUE)
  if (length(mult_cols) > 0) {
    for (col in mult_cols) {
      level_name <- sub("^effect_multiplier_", "", col)
      posterior_name <- paste0("effect_multiplier[", level_name, "]")
      param_map[[posterior_name]] <- col
    }
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
      min_ess_bulk    = post_row$min_ess_bulk,
      # Per-fit sampler diagnostics. Same value for every row of a given
      # fit (one fit = one rep); aggregation below pulls a single value
      # per rep before computing rep-level stats.
      n_divergent     = post_row$n_divergent     %||% NA_integer_,
      n_max_treedepth = post_row$n_max_treedepth %||% NA_integer_,
      min_ebfmi       = post_row$min_ebfmi       %||% NA_real_
    )
  })

  results
}

# Local %||% fallback so this file is independent of rlang attach state.
`%||%` <- function(x, y) if (is.null(x)) y else x

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
      # Monte Carlo standard errors, so the figure can carry uncertainty
      # without going back to rep-level data. These are the MC error of the
      # summary across replicates (sd / sqrt(n)), NOT posterior uncertainty.
      se_bias            = sd(bias, na.rm = TRUE) / sqrt(sum(!is.na(bias))),
      se_ci_width        = sd(ci_width, na.rm = TRUE) /
                             sqrt(sum(!is.na(ci_width))),
      se_posterior_sd    = sd(posterior_sd, na.rm = TRUE) /
                             sqrt(sum(!is.na(posterior_sd))),
      # Convergence diagnostics
      max_rhat_worst     = max(max_rhat, na.rm = TRUE),
      min_ess_worst      = min(min_ess_bulk, na.rm = TRUE),
      # Fraction with Rhat > 1.05 (problematic fits)
      pct_rhat_bad       = mean(max_rhat > 1.05, na.rm = TRUE),
      # Sampler diagnostics (n_divergent etc. are constant within a rep
      # so the summary across rows-of-this-rep equals the value itself;
      # aggregating across reps gives a per-scenario picture).
      mean_n_divergent           = mean(n_divergent, na.rm = TRUE),
      max_n_divergent            = max(n_divergent,  na.rm = TRUE),
      pct_reps_with_divergence   = mean(n_divergent > 0,  na.rm = TRUE),
      pct_reps_many_divergence   = mean(n_divergent > 10, na.rm = TRUE),
      mean_n_max_treedepth       = mean(n_max_treedepth, na.rm = TRUE),
      min_ebfmi_worst            = min(min_ebfmi, na.rm = TRUE),
      .groups = "drop"
    )
}

# ===========================================================================
# Paired RMSE contrast between two model arms
# ===========================================================================

#' Paired full/naive RMSE ratio with a Monte Carlo interval
#'
#' Both comparison arms of a scenario are fitted to the SAME simulated dataset
#' within a replicate (see run_one_rep()), so their errors are correlated.
#' Taking the ratio of two separately-aggregated RMSEs throws that pairing
#' away and inflates the Monte Carlo error substantially; resampling REPLICATES
#' (not arms) keeps it. At the low replication counts the figure sweeps use,
#' this is the difference between a readable curve and noise.
#'
#' @param rep_results Per-replication results with model_label, bias, and the
#'   tar_map_rep replicate keys (tar_batch / tar_rep).
#' @param num,den     Model labels for the numerator and denominator.
#' @param param       Parameter to contrast.
#' @param n_boot      Paired bootstrap resamples.
#' @param seed        Fixed for reproducibility across pipeline runs.
#' @return Tibble: scenario_id, n_reps, ratio, lo, hi.
paired_ratio <- function(rep_results, column = "bias",
                         stat = c("rmse", "mean"),
                         num = "full", den = "naive",
                         param = "treatment_effect_mean",
                         n_boot = 2000L, seed = 1L) {
  stat <- match.arg(stat)
  keys <- intersect(c("tar_batch", "tar_rep"), names(rep_results))
  if (!length(keys)) {
    stop("rep_results has no tar_batch/tar_rep replicate key")
  }

  wide <- rep_results |>
    dplyr::filter(parameter == param, model_label %in% c(num, den)) |>
    dplyr::select(dplyr::all_of(c("scenario_id", keys, "model_label", column))) |>
    tidyr::pivot_wider(names_from = model_label,
                       values_from = dplyr::all_of(column))

  if (!all(c(num, den) %in% names(wide))) return(tibble::tibble())

  ratio_fn <- switch(stat,
    rmse = function(a, b) sqrt(mean(a^2) / mean(b^2)),
    mean = function(a, b) mean(a) / mean(b))

  withr::with_seed(seed, {
    wide |>
      dplyr::group_by(scenario_id) |>
      dplyr::group_modify(function(g, ...) {
        a <- g[[num]]; b <- g[[den]]
        ok <- is.finite(a) & is.finite(b)
        a <- a[ok]; b <- b[ok]
        if (length(a) < 2L) {
          return(tibble::tibble(n_reps = length(a), ratio = NA_real_,
                                lo = NA_real_, hi = NA_real_))
        }
        # Resample REPLICATES, carrying both arms together -- this is what
        # preserves the pairing.
        boot <- vapply(seq_len(n_boot), function(i) {
          idx <- sample.int(length(a), replace = TRUE)
          ratio_fn(a[idx], b[idx])
        }, numeric(1))
        tibble::tibble(
          n_reps = length(a),
          ratio  = ratio_fn(a, b),
          lo     = unname(quantile(boot, 0.05, na.rm = TRUE)),
          hi     = unname(quantile(boot, 0.95, na.rm = TRUE))
        )
      }) |>
      dplyr::ungroup()
  })
}
