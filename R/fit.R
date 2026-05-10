# Model fitting wrappers for the metadid validation study.
#
# fit_scenario()    — fit a single model configuration to simulated data
# fit_comparison()  — fit multiple model configurations to the same data (for B-category)

library(metadid)

# MCMC configuration (shared across all fits)
MCMC_OPTS <- list(
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 1000L,
  refresh         = 0
)

# ===========================================================================
# Single-model fit
# ===========================================================================

#' Fit a model to simulated data according to scenario config
#'
#' @param sim_result Output from simulate_scenario() (list with $data, $true_params)
#' @param fit_config A fit config list (from scenario$fit or a compare element)
#' @return A tibble of posterior summaries for key parameters
fit_scenario <- function(sim_result, fit_config) {

  data <- sim_result$data

  # Build argument list
  args <- c(
    data,  # unpacks to summary_data = ... and/or individual_data = ...
    list(
      normalise_by_baseline = fit_config$normalise_by_baseline,
      robust_heterogeneity  = fit_config$robust_heterogeneity,
      design_effects        = fit_config$design_effects,
      hierarchical_rho      = fit_config$hierarchical_rho
    ),
    MCMC_OPTS
  )

  # Add covariates if specified
  if (!is.null(fit_config$covariates)) {
    args$covariates <- fit_config$covariates
  }

  # Correlated effects (only for meta_did_general)
  if (isTRUE(fit_config$correlated_effects)) {
    args$correlated_effects <- TRUE
  }

  # meta_did_general-specific options
  if (identical(fit_config$fn, "meta_did_general")) {
    args$time_trend         <- fit_config$time_trend
    args$baseline_imbalance <- fit_config$baseline_imbalance
    args$pp_likelihood      <- fit_config$pp_likelihood
  }

  # Call the fitting function
  fit_fn <- switch(
    fit_config$fn,
    meta_did         = metadid::meta_did,
    meta_did_general = metadid::meta_did_general,
    stop("Unknown fit function: ", fit_config$fn)
  )
  fit <- do.call(fit_fn, args)

  # Extract posterior summaries
  extract_posteriors(fit)
}

# ===========================================================================
# Posterior extraction
# ===========================================================================

#' Extract posterior summaries from a meta_did_fit object
#'
#' Returns a tibble with columns: parameter, mean, sd, lo, hi
#' where lo/hi are the 5th/95th percentiles (90% CI).
#'
#' Supplements the S3 summary() with additional Stan parameters
#' (time_trend_mean, time_trend_sd, baseline_control_mean) that
#' the summary method does not expose.
extract_posteriors <- function(fit) {
  s <- summary(fit, prob = 0.9)

  # Additional parameters to extract directly from Stan draws
  extra_params <- intersect(
    c("time_trend_mean", "time_trend_sd", "baseline_control_mean"),
    fit$fit$metadata()$stan_variables
  )

  if (length(extra_params) > 0 && fit$method == "sample") {
    lo_q <- 0.05
    hi_q <- 0.95
    extra_rows <- purrr::map_dfr(extra_params, function(p) {
      d <- fit$fit$draws(p, format = "matrix")
      tibble::tibble(
        parameter = p,
        mean      = mean(d),
        sd        = sd(d),
        lo        = unname(quantile(d, lo_q)),
        hi        = unname(quantile(d, hi_q))
      )
    })
    s <- dplyr::bind_rows(s, extra_rows)
  }

  # Convergence diagnostics for key population parameters
  key_params <- intersect(
    c("treatment_effect_mean", "treatment_effect_sd",
      "time_trend_mean", "time_trend_sd",
      "baseline_control_mean"),
    fit$fit$metadata()$stan_variables
  )

  if (length(key_params) > 0 && fit$method == "sample") {
    diag_df <- fit$fit$summary(key_params)
    max_rhat <- max(diag_df$rhat, na.rm = TRUE)
    min_ess  <- min(diag_df$ess_bulk, na.rm = TRUE)
  } else {
    max_rhat <- NA_real_
    min_ess  <- NA_real_
  }

  s |>
    dplyr::mutate(
      max_rhat = max_rhat,
      min_ess_bulk = min_ess
    )
}

# ===========================================================================
# Run one replication (simulate + fit + assess)
# ===========================================================================

#' Run a single replication: simulate, fit, assess
#'
#' For standard scenarios (no compare), simulates data, fits one model,
#' and returns assessment results.
#'
#' For comparative scenarios, simulates once and fits all model configs,
#' returning assessment results with a model_label column.
run_one_rep <- function(scenario_id, rep_seed) {
  config <- SCENARIO_CONFIGS[[scenario_id]]
  sim    <- simulate_scenario(scenario_id, rep_seed)

  if (!is.null(config$compare)) {
    # Comparative study: fit each config
    results <- lapply(config$compare, function(cmp) {
      # Merge compare overrides into base fit config
      fit_config <- modifyList(config$fit, cmp)

      # Handle data_format override (may change summary vs individual)
      if (!is.null(cmp$data_format) && cmp$data_format != config$fit$data_format) {
        # Re-assemble data with the overridden format
        fit_config_for_data <- modifyList(config$fit, cmp)
        sim_data <- reassemble_data(sim, config$dgp, fit_config_for_data)
      } else {
        sim_data <- sim
      }

      posteriors <- fit_scenario(sim_data, fit_config)
      assessment <- assess_one(posteriors, sim$true_params, fit_config)
      assessment$model_label <- cmp$label
      assessment
    })
    result <- dplyr::bind_rows(results)
  } else {
    # Standard: single fit
    posteriors <- fit_scenario(sim, config$fit)
    result     <- assess_one(posteriors, sim$true_params, config$fit)
    result$model_label <- "default"
  }

  result$scenario_id <- scenario_id
  result
}

# Helper: re-assemble data when a comparison needs a different format
reassemble_data <- function(sim, dgp, fit_config) {
  # For individual vs summary comparisons, we need to re-simulate

  # since we can't easily reconstruct individual data from summary.
  # Instead, we re-run simulation with the same seed — but this is tricky.
  # A simpler approach: for A12, the simulate_scenario already returns
  # the raw simulation; we just re-assemble.
  #
  # For now, this only handles the provide_rho override case.
  sim_copy <- sim
  if (!is.null(sim_copy$data$summary_data) && !fit_config$provide_rho) {
    sim_copy$data$summary_data$rho <- NA_real_
  }
  sim_copy
}
