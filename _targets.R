# _targets.R — metadid simulation-based validation pipeline
#
# Run the full pipeline:
#   targets::tar_make()
#
# Run individual categories:
#   targets::tar_make(names = starts_with("A_"))
#   targets::tar_make(names = starts_with("F_"))
#   targets::tar_make(names = starts_with("B_"))
#   targets::tar_make(names = starts_with("C_"))
#   targets::tar_make(names = starts_with("D_"))
#   targets::tar_make(names = starts_with("E_"))
#
# Render the validation report:
#   targets::tar_make(names = "report")

library(targets)
library(tarchetypes)

# ---------------------------------------------------------------------------
# Global configuration
# ---------------------------------------------------------------------------

N_REPS <- 25L

# ---------------------------------------------------------------------------
# Source R files
# ---------------------------------------------------------------------------

tar_source("R/scenarios.R")
tar_source("R/simulate.R")
tar_source("R/fit.R")
tar_source("R/assess.R")
tar_source("R/plots.R")

# ---------------------------------------------------------------------------
# Pipeline definition
# ---------------------------------------------------------------------------

list(
  # ---- Category A: Calibration studies ----
  tarchetypes::tar_map_rep(
    name    = A_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("A")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(A_agg, aggregate_scenario(A_rep)),

  # ---- Category F: Large-N bias probes ----
  tarchetypes::tar_map_rep(
    name    = F_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("F")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(F_agg, aggregate_scenario(F_rep)),

  # ---- Category B: Comparative studies ----
  tarchetypes::tar_map_rep(
    name    = B_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("B")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(B_agg, aggregate_scenario(B_rep)),

  # ---- Category C: Outlier and heavy-tailed ----
  tarchetypes::tar_map_rep(
    name    = C_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("C")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(C_agg, aggregate_scenario(C_rep)),

  # ---- Category D: Assumption violations ----
  tarchetypes::tar_map_rep(
    name    = D_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("D")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(D_agg, aggregate_scenario(D_rep)),

  # ---- Category E: Edge cases ----
  tarchetypes::tar_map_rep(
    name    = E_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get()),
    values  = tibble::tibble(scenario_id = scenario_ids("E")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(E_agg, aggregate_scenario(E_rep)),

  # ---- Combined results ----
  tar_target(
    all_agg,
    dplyr::bind_rows(A_agg, F_agg, B_agg, C_agg, D_agg, E_agg)
  ),

  tar_target(
    all_rep,
    dplyr::bind_rows(A_rep, F_rep, B_rep, C_rep, D_rep, E_rep)
  ),

  # ---- Machine-readable exports ----
  tar_target(
    export_agg_csv,
    {
      dir.create("output", showWarnings = FALSE)
      readr::write_csv(all_agg, "output/aggregated_results.csv")
      "output/aggregated_results.csv"
    },
    format = "file"
  ),

  tar_target(
    export_rep_rds,
    {
      dir.create("output", showWarnings = FALSE)
      saveRDS(all_rep, "output/replication_results.rds")
      "output/replication_results.rds"
    },
    format = "file"
  ),

  # ---- Validation report ----
  tarchetypes::tar_quarto(
    report,
    path = "reports/validation-report.qmd"
  )
)
