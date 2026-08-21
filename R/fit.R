# Model fitting wrappers for the metadid validation study.
#
# fit_scenario()    — fit a single model configuration to simulated data
# fit_comparison()  — fit multiple model configurations to the same data (for B-category)

library(metadid)

# MCMC configuration (shared across all fits).
#
# chains = 2, parallel_chains = 2 means each fit runs its 2 chains
# concurrently on 2 cores in time T (single-chain wall time). The targets
# controller in _targets.R sets workers = 2, so two fits run side by side
# on the 4-core GitHub Actions runner — 2 fits × 2 chains = 4 chains in
# flight, completing a pair of fits in T versus the previous 2T for two
# sequential 4-chain fits. Net: roughly 2x throughput per scenario.
#
# The cost is statistical: 2 chains is the minimum for Rhat-style
# convergence diagnostics, and gives weaker signal on between-chain
# disagreement than the conventional 4 chains. If divergences or Rhat
# become a problem, bump chains back up.
MCMC_OPTS <- list(
  chains          = 2L,
  parallel_chains = 2L,
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
      normalise_by_baseline = fit_config$normalise,
      robust_heterogeneity  = fit_config$robust_heterogeneity,
      design_effects        = fit_config$design_effects,
      hierarchical_rho      = fit_config$hierarchical_rho
    ),
    MCMC_OPTS
  )

  # Custom priors, when a scenario (or a compare arm) supplies them. Used by
  # the M-category `homog` arm, which pins time_trend_sd at ~0 to emulate a
  # homogeneous-trend model — metadid's `time_trend` switch only offers
  # "pooled" and "fixed_zero", so there is no built-in sigma_beta = 0 option.
  if (!is.null(fit_config$priors)) {
    args$priors <- fit_config$priors
  }

  # Add covariates if specified
  if (!is.null(fit_config$covariates)) {
    args$covariates <- fit_config$covariates
  }

  # Add multiplicative covariate if specified (multiplicative-covars feature).
  # Either a single column name (character) or a one-sided formula naming up
  # to two columns (~ a + b, product structure). Studies at non-reference
  # levels have their population-mean linear predictor multiplied by the
  # estimated effect_multiplier[<level>] (one covariate) or
  # effect_multiplier[<column>:<level>] factors (two covariates). Used by
  # the I-category scenarios.
  if (!is.null(fit_config$multiplicative_covariate)) {
    args$multiplicative_covariate <- fit_config$multiplicative_covariate
  }

  # Correlated effects (only for meta_did_general)
  if (isTRUE(fit_config$correlated_effects)) {
    args$correlated_effects <- TRUE
  }

  # Baseline-imbalance options. Both meta_did() and meta_did_general() accept
  # these, so they are passed regardless of which entry point an arm uses --
  # that is what lets a "full" and a "naive" arm differ only in gamma handling.
  args$baseline_imbalance   <- fit_config$baseline_imbalance   %||% "by_randomisation"
  args$mu_gamma             <- fit_config$mu_gamma             %||% "zero"
  args$kappa                <- fit_config$kappa                %||% 0.5
  args$cluster_deff_default <- fit_config$cluster_deff_default %||% 2
  args$allow_unidentified_kappa <- fit_config$allow_unidentified_kappa %||% FALSE

  # meta_did_general-specific options
  if (identical(fit_config$fn, "meta_did_general")) {
    args$time_trend    <- fit_config$time_trend
    args$pp_likelihood <- fit_config$pp_likelihood
  }

  # Thousands of fits run headless; the per-fit advisory about a missing
  # randomisation column would drown the logs. The scenarios set the column
  # deliberately, so its absence is a choice here rather than an oversight.
  old_quiet <- getOption("metadid.quiet")
  options(metadid.quiet = TRUE)
  on.exit(options(metadid.quiet = old_quiet), add = TRUE)

  # Call the fitting function
  fit_fn <- switch(
    fit_config$fn,
    meta_did         = metadid::meta_did,
    meta_did_general = metadid::meta_did_general,
    stop("Unknown fit function: ", fit_config$fn)
  )

  # A few percent of fits die reading their own CmdStan output -- "File does
  # not exist: .../meta_analysis_master-<stamp>-<chain>-<hash>.csv", or a
  # truncated CSV surfacing as "length of 'dimnames' [2] not equal to array
  # extent". Measured at 16/400 (4%) locally, running SERIALLY on an idle
  # machine, so it is not contention between crew workers. In CI a single such
  # branch kills the whole multi-hour job.
  #
  # The retry deliberately re-runs with IDENTICAL arguments, including the
  # seed. That matters for the statistics: a successful retry is the same
  # posterior draw the first attempt was computing, so nothing is selected on.
  # Retrying with a fresh seed would instead condition on success and could
  # quietly bias results toward well-behaved datasets. A failure that is really
  # about the data or the model reproduces under the same seed and still
  # surfaces, which is what we want.
  .attempt <- function(k) {
    fit <- do.call(fit_fn, args)
    extract_posteriors(fit)
  }
  # 5, not 3. The camera-ready run is ~6,280 fits (157 per rep x 40 reps), so
  # at a 4% per-fit flake rate P(some fit exhausts its attempts, killing the
  # job) is 1.00 at 2 attempts, 0.33 at 3, 0.016 at 4 and 0.001 at 5. Retries
  # only fire on the 4%, so the expected extra compute is ~4%.
  n_try <- as.integer(Sys.getenv("FIT_MAX_ATTEMPTS", "5"))
  for (k in seq_len(n_try)) {
    out <- try(.attempt(k), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    msg <- conditionMessage(attr(out, "condition"))
    if (k == n_try) {
      stop(sprintf("fit failed after %d attempts with identical inputs: %s",
                   n_try, msg), call. = FALSE)
    }
    message(sprintf("fit attempt %d/%d failed (%s); retrying with the same seed",
                    k, n_try, msg))
  }
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
    c("time_trend_mean", "time_trend_sd", "baseline_control_mean",
      # Baseline-imbalance block. baseline_difference_mean and kappa are
      # transformed parameters that exist whether or not they are sampled, so
      # they always resolve; a pinned one has zero posterior SD and assess_one()
      # drops it rather than crediting a constant with perfect recovery.
      "baseline_difference_mean", "baseline_difference_sd", "kappa"),
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

  # Sampler diagnostics aggregated across chains for the whole fit.
  # diagnostic_summary() returns per-chain vectors; we sum / max as
  # appropriate. Suppressing the warning chatter — the same information
  # comes back in the values.
  if (fit$method == "sample") {
    diag <- suppressWarnings(fit$fit$diagnostic_summary(quiet = TRUE))
    n_divergent      <- sum(diag$num_divergent      %||% 0L)
    n_max_treedepth  <- sum(diag$num_max_treedepth  %||% 0L)
    min_ebfmi        <- min(diag$ebfmi              %||% NA_real_, na.rm = TRUE)
  } else {
    n_divergent     <- NA_integer_
    n_max_treedepth <- NA_integer_
    min_ebfmi       <- NA_real_
  }

  s |>
    dplyr::mutate(
      max_rhat         = max_rhat,
      min_ess_bulk     = min_ess,
      n_divergent      = n_divergent,
      n_max_treedepth  = n_max_treedepth,
      min_ebfmi        = min_ebfmi
    )
}

# Local %||% fallback in case rlang isn't attached at use time.
`%||%` <- function(x, y) if (is.null(x)) y else x

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
run_one_rep <- function(scenario_id, config, rep_seed, pkg = NULL) {
  sim <- simulate_scenario(scenario_id, rep_seed, config)

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
  if (!is.null(sim_copy$data$summary_data) &&
      "rho" %in% names(sim_copy$data$summary_data)) {
    sim_copy$data$summary_data$rho <-
      apply_provide_rho(sim_copy$data$summary_data$rho, fit_config$provide_rho)
  }
  sim_copy
}
