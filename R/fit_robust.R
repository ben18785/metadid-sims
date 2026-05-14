# fit_robust.R — robust model comparator fits for A and F categories.
# Kept separate from fit.R so that changes here do not invalidate A or F targets.

library(metadid)

#' Run one replication for a scenario with robust_heterogeneity = TRUE
#'
#' Simulates the same data as run_one_rep() would for the same scenario_id
#' and rep_seed, then fits the robust model.
#'
#' @param scenario_id String identifying the scenario in SCENARIO_CONFIGS
#' @param rep_seed    Integer seed for this replication
#' @param pkg         Ignored; used to create a dependency on metadid_src
run_robust_rep <- function(scenario_id, rep_seed, pkg = NULL) {

  config        <- SCENARIO_CONFIGS[[scenario_id]]
  sim           <- simulate_scenario(scenario_id, rep_seed)

  robust_config                        <- config
  robust_config$fit$robust_heterogeneity <- TRUE

  posts  <- fit_scenario(sim, robust_config$fit)
  result <- assess_one(posts, sim$true_params, robust_config$fit)

  result$model_label <- "robust"
  result$scenario_id <- scenario_id
  result
}
