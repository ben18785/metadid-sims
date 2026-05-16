# fit_g.R — fitting wrapper for Category G bias-source investigation scenarios.
# Kept separate from fit.R so that changes here do not invalidate A–F targets.

library(metadid)

#' Run one replication for a G scenario, supporting custom priors
#'
#' Follows the same structure as run_one_rep() but passes fit_config$priors
#' to the fitting function when present.
#'
#' @param scenario_id String identifying the scenario in SCENARIO_CONFIGS
#' @param rep_seed    Integer seed for this replication
#' @param pkg         Ignored; used to create a dependency on metadid_src
run_g_rep <- function(scenario_id, config, rep_seed, pkg = NULL) {

  sim <- simulate_scenario(scenario_id, rep_seed, config)

  # Build fit args, identical to fit_scenario() but also passing priors
  fit_config <- config$fit
  args <- c(
    sim$data,
    list(
      normalise_by_baseline = fit_config$normalise_by_baseline,
      robust_heterogeneity  = fit_config$robust_heterogeneity,
      design_effects        = fit_config$design_effects,
      hierarchical_rho      = fit_config$hierarchical_rho
    ),
    MCMC_OPTS
  )

  if (!is.null(fit_config$priors)) {
    args$priors <- fit_config$priors
  }

  fit_fn <- switch(
    fit_config$fn,
    meta_did         = meta_did,
    meta_did_general = meta_did_general,
    stop("Unknown fit function: ", fit_config$fn)
  )

  fit    <- do.call(fit_fn, args)
  posts  <- extract_posteriors(fit)

  # Strip priors before assess_one so it doesn't encounter the complex list
  fit_config_clean        <- fit_config
  fit_config_clean$priors <- NULL

  result             <- assess_one(posts, sim$true_params, fit_config_clean)
  result$model_label <- "default"
  result$scenario_id <- scenario_id
  result
}
