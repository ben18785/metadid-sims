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
tar_source("R/fit_g.R")
tar_source("R/fit_robust.R")
tar_source("R/assess.R")
tar_source("R/plots.R")

# ---------------------------------------------------------------------------
# Pipeline definition
# ---------------------------------------------------------------------------

list(
  # ---- Track metadid source: invalidates all downstream targets on reinstall ----
  tar_target(
    metadid_src,
    c(list.files("../metadid/R", full.names = TRUE), "../metadid/DESCRIPTION"),
    format = "file"
  ),

  # ---- Category A: Calibration studies ----
  tarchetypes::tar_map_rep(
    name    = A_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("A")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(A_agg, aggregate_scenario(A_rep)),

  # ---- Category A: Robust comparator (new fits, does not invalidate A_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust_heterogeneity)
  tarchetypes::tar_map_rep(
    name    = A_rep_robust,
    command = run_robust_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = Filter(
      function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
      scenario_ids("A")
    )),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(
    A_combined,
    dplyr::bind_rows(
      dplyr::mutate(A_rep,        model_label = "normal"),
      dplyr::mutate(A_rep_robust, model_label = "robust")
    )
  ),
  tar_target(A_agg_combined, aggregate_scenario(A_combined)),

  # ---- Category F: Large-N bias probes ----
  tarchetypes::tar_map_rep(
    name    = F_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("F")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(F_agg, aggregate_scenario(F_rep)),

  # ---- Category F: Robust comparator (new fits, does not invalidate F_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust_heterogeneity)
  tarchetypes::tar_map_rep(
    name    = F_rep_robust,
    command = run_robust_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = Filter(
      function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
      scenario_ids("F")
    )),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(
    F_combined,
    dplyr::bind_rows(
      dplyr::mutate(F_rep,        model_label = "normal"),
      dplyr::mutate(F_rep_robust, model_label = "robust")
    )
  ),
  tar_target(F_agg_combined, aggregate_scenario(F_combined)),

  # ---- Category B: Comparative studies ----
  tarchetypes::tar_map_rep(
    name    = B_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("B")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(B_agg, aggregate_scenario(B_rep)),

  # ---- Category C: Outlier and heavy-tailed ----
  tarchetypes::tar_map_rep(
    name    = C_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("C")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(C_agg, aggregate_scenario(C_rep)),

  # ---- Category D: Assumption violations ----
  tarchetypes::tar_map_rep(
    name    = D_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("D")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(D_agg, aggregate_scenario(D_rep)),

  # ---- Category E: Edge cases ----
  tarchetypes::tar_map_rep(
    name    = E_rep,
    command = run_one_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("E")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(E_agg, aggregate_scenario(E_rep)),

  # ---- Scenario lookup table ----
  tar_target(scenario_lookup_tbl, scenario_lookup()),

  # ---- Category G: Bias source investigation ----
  # Uses run_g_rep() (from fit_g.R) rather than run_one_rep() so that changes
  # to the priors-aware fitting code do not invalidate A-F targets.
  tarchetypes::tar_map_rep(
    name    = G_rep,
    command = run_g_rep(scenario_id, targets::tar_seed_get(), metadid_src),
    values  = tibble::tibble(scenario_id = scenario_ids("G")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1L
  ),
  tar_target(G_agg, aggregate_scenario(G_rep)),

  # ---- Diagnostic plots (% bias, coverage, Rhat) per category ----
  tar_target(diag_plot_A, plot_diagnostics(A_rep, A_agg, scenario_lookup_tbl, "A")),
  tar_target(diag_plot_G, plot_diagnostics(G_rep, G_agg, scenario_lookup_tbl, "G")),
  tar_target(diag_plot_F, plot_diagnostics(F_rep, F_agg, scenario_lookup_tbl, "F")),
  tar_target(diag_plot_B, plot_diagnostics(B_rep, B_agg, scenario_lookup_tbl, "B")),
  tar_target(diag_plot_C, plot_diagnostics(C_rep, C_agg, scenario_lookup_tbl, "C")),
  tar_target(diag_plot_D, plot_diagnostics(D_rep, D_agg, scenario_lookup_tbl, "D")),
  tar_target(diag_plot_E, plot_diagnostics(E_rep, E_agg, scenario_lookup_tbl, "E")),

  # ---- Combined results ----
  tar_target(
    all_agg,
    dplyr::bind_rows(A_agg, F_agg, B_agg, C_agg, D_agg, E_agg, G_agg)
  ),

  tar_target(
    all_rep,
    dplyr::bind_rows(A_rep, F_rep, B_rep, C_rep, D_rep, E_rep, G_rep)
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
    archive_results,
    {
      install_datetime <- file.info(
        system.file("DESCRIPTION", package = "metadid")
      )$mtime

      stamped <- all_agg |>
        dplyr::mutate(
          metadid_install_datetime = install_datetime,
          run_date                 = Sys.time()
        )

      archive_path <- "output/archive.csv"
      if (file.exists(archive_path)) {
        existing <- readr::read_csv(archive_path, show_col_types = FALSE)
        # Only append if this install datetime isn't already recorded
        if (!any(existing$metadid_install_datetime == install_datetime)) {
          stamped <- dplyr::bind_rows(existing, stamped)
        } else {
          stamped <- existing
        }
      }

      readr::write_csv(stamped, archive_path)
      archive_path
    },
    format = "file",
    cue = tar_cue(mode = "always")
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
