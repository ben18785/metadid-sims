# fit_naive.R — naive model comparator fits for Category A.
# Kept separate from fit.R so that changes here do not invalidate A or other targets.
#
# "Naive" means meta_did_general with time_trend = "fixed_zero" and
# baseline_imbalance = "fixed_zero": the simplest possible pooled model that
# ignores time trends and baseline differences between arms.

library(metadid)

#' Run one replication for an A scenario with the naive model
#'
#' Simulates the same data as run_one_rep() for the same scenario_id and
#' rep_seed, then fits meta_did_general with fixed-zero trend and imbalance.
#' normalise_by_baseline, design_effects, and hierarchical_rho are preserved
#' from the base scenario config so like-for-like comparisons are valid.
#'
#' @param scenario_id String identifying the scenario in SCENARIO_CONFIGS
#' @param rep_seed    Integer seed for this replication
#' @param pkg         Ignored; used to create a dependency on metadid_src
run_naive_rep <- function(scenario_id, config, rep_seed, pkg = NULL) {

  sim <- simulate_scenario(scenario_id, rep_seed, config)

  # Build naive config: meta_did_general, zero trend, zero imbalance
  naive_config <- list(
    fn                    = "meta_did_general",
    normalise_by_baseline = config$fit$normalise_by_baseline,
    robust_heterogeneity  = FALSE,
    design_effects        = config$fit$design_effects,
    hierarchical_rho      = config$fit$hierarchical_rho,
    correlated_effects    = FALSE,
    time_trend            = "fixed_zero",
    baseline_imbalance    = "fixed_zero",
    data_format           = "summary"
  )

  # Naive model always uses summary data
  naive_sim      <- sim
  naive_sim$data <- list(summary_data = sim$data$summary_data)

  posts  <- fit_scenario(naive_sim, naive_config)
  result <- assess_one(posts, sim$true_params, naive_config)

  result$model_label <- "naive"
  result$scenario_id <- scenario_id
  result
}
