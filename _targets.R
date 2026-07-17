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
#   targets::tar_make(names = starts_with("H_"))
#   targets::tar_make(names = starts_with("G_"))
#
# Render the validation report:
#   targets::tar_make(names = "report")

library(targets)
library(tarchetypes)

# ---------------------------------------------------------------------------
# Parallel execution
# ---------------------------------------------------------------------------
#
# Run two targets concurrently within each pipeline invocation. Combined
# with parallel_chains = 2 in R/fit.R, this saturates the 4-core GitHub
# Actions runner: 2 concurrent fits × 2 chains per fit = 4 chains in
# flight. Reduces wall time within a scenario when there are many
# replications/fit-variants (e.g. A4 = 3 fit-variants × 25 reps = 75
# sequential fits previously).
#
# crew_controller_local launches local R subprocesses; no SSH/cluster
# setup required.
tar_option_set(
  controller = crew::crew_controller_local(workers = 2L)
)

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
tar_source("R/fit_naive.R")
tar_source("R/assess.R")
tar_source("R/plots.R")

# ---------------------------------------------------------------------------
# Helper: build a standard values tibble for tar_map_rep
#
# Embeds each scenario's config directly in the values so that targets tracks
# per-scenario dependencies rather than the monolithic SCENARIO_CONFIGS list.
# Changing one scenario's DGP or fit config now only invalidates that scenario's
# targets; all others remain cached.
# ---------------------------------------------------------------------------

scenario_values <- function(ids) {
  tibble::tibble(
    scenario_id = ids,
    config      = lapply(ids, function(s) SCENARIO_CONFIGS[[s]])
  )
}

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
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("A")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(A_agg, aggregate_scenario(A_rep)),

  # ---- Category A: Naive comparator (does not invalidate A_rep) ----
  # Excludes scenarios with data_format = "individual" (no summary_data available)
  tarchetypes::tar_map_rep(
    name    = A_rep_naive,
    command = run_naive_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$data_format == "individual"),
        scenario_ids("A")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),

  # ---- Category A: Robust comparator (does not invalidate A_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust)
  tarchetypes::tar_map_rep(
    name    = A_rep_robust,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
        scenario_ids("A")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(
    A_combined,
    dplyr::bind_rows(
      dplyr::mutate(A_rep,        model_label = "normal"),
      dplyr::mutate(A_rep_robust, model_label = "robust"),
      A_rep_naive  # model_label = "naive" set in run_naive_rep()
    )
  ),
  tar_target(A_agg_combined, aggregate_scenario(A_combined)),

  # ---- Category F: Large-N bias probes ----
  tarchetypes::tar_map_rep(
    name    = F_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("F")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(F_agg, aggregate_scenario(F_rep)),

  # ---- Category F: Robust comparator (does not invalidate F_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust)
  tarchetypes::tar_map_rep(
    name    = F_rep_robust,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
        scenario_ids("F")
      )
      scenario_values(ids)
    },
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
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("B")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(B_agg, aggregate_scenario(B_rep)),

  # ---- Category C: Outlier and heavy-tailed ----
  tarchetypes::tar_map_rep(
    name    = C_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("C")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(C_agg, aggregate_scenario(C_rep)),

  # ---- Category D: Assumption violations ----
  tarchetypes::tar_map_rep(
    name    = D_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("D")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(D_agg, aggregate_scenario(D_rep)),

  # ---- Category D: Naive comparator (does not invalidate D_rep) ----
  # Excludes scenarios with data_format = "individual" (no summary_data available)
  tarchetypes::tar_map_rep(
    name    = D_rep_naive,
    command = run_naive_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$data_format == "individual"),
        scenario_ids("D")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),

  # ---- Category D: Robust comparator (does not invalidate D_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust)
  tarchetypes::tar_map_rep(
    name    = D_rep_robust,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
        scenario_ids("D")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(
    D_combined,
    dplyr::bind_rows(
      dplyr::mutate(D_rep,        model_label = "normal"),
      dplyr::mutate(D_rep_robust, model_label = "robust"),
      D_rep_naive  # model_label = "naive" set in run_naive_rep()
    )
  ),
  tar_target(D_agg_combined, aggregate_scenario(D_combined)),

  # ---- Category E: Edge cases ----
  tarchetypes::tar_map_rep(
    name    = E_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("E")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(E_agg, aggregate_scenario(E_rep)),

  # ---- Category E: Naive comparator (does not invalidate E_rep) ----
  # Excludes scenarios with data_format = "individual" (no summary_data available)
  tarchetypes::tar_map_rep(
    name    = E_rep_naive,
    command = run_naive_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$data_format == "individual"),
        scenario_ids("E")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),

  # ---- Category E: Robust comparator (does not invalidate E_rep) ----
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust)
  tarchetypes::tar_map_rep(
    name    = E_rep_robust,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = {
      ids <- Filter(
        function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects),
        scenario_ids("E")
      )
      scenario_values(ids)
    },
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(
    E_combined,
    dplyr::bind_rows(
      dplyr::mutate(E_rep,        model_label = "normal"),
      dplyr::mutate(E_rep_robust, model_label = "robust"),
      E_rep_naive  # model_label = "naive" set in run_naive_rep()
    )
  ),
  tar_target(E_agg_combined, aggregate_scenario(E_combined)),

  # ---- Category H: Time trend distributional misspecification ----
  tarchetypes::tar_map_rep(
    name    = H_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("H")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(H_agg, aggregate_scenario(H_rep)),

  # ---- Category I: Multiplicative covariate scenarios ----
  # Each scenario handles its own modelled/raw and with/without comparators
  # via the per-scenario `compare` block in scenarios.R, so no separate
  # *_naive / *_robust comparator branches are needed.
  tarchetypes::tar_map_rep(
    name    = I_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("I")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tar_target(I_agg, aggregate_scenario(I_rep)),

  # ---- Scenario lookup table ----
  tar_target(scenario_lookup_tbl, scenario_lookup()),

  # ---- Category G: Bias source investigation ----
  # Uses run_g_rep() (from fit_g.R) so that changes to priors-aware fitting
  # code do not invalidate A-F targets.
  tarchetypes::tar_map_rep(
    name    = G_rep,
    command = run_g_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("G")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS * 2,
    reps    = 1L
  ),
  tar_target(G_agg, aggregate_scenario(G_rep)),

  # ---- Diagnostic plots (% bias, coverage, Rhat) per category ----
  tar_target(diag_plot_A, plot_diagnostics(A_rep,  A_agg,  scenario_lookup_tbl, "A")),
  tar_target(diag_plot_F, plot_diagnostics(F_rep,  F_agg,  scenario_lookup_tbl, "F")),
  tar_target(diag_plot_B, plot_diagnostics(B_rep,  B_agg,  scenario_lookup_tbl, "B")),
  tar_target(diag_plot_C, plot_diagnostics(C_rep,  C_agg,  scenario_lookup_tbl, "C")),
  tar_target(diag_plot_D, plot_diagnostics(D_rep,  D_agg,  scenario_lookup_tbl, "D")),
  tar_target(diag_plot_E, plot_diagnostics(E_rep,  E_agg,  scenario_lookup_tbl, "E")),
  tar_target(diag_plot_G, plot_diagnostics(G_rep,  G_agg,  scenario_lookup_tbl, "G")),
  tar_target(diag_plot_H, plot_diagnostics(H_rep,  H_agg,  scenario_lookup_tbl, "H")),
  tar_target(diag_plot_I, plot_diagnostics(I_rep,  I_agg,  scenario_lookup_tbl, "I")),

  # ---- Combined results ----
  tar_target(
    all_agg,
    dplyr::bind_rows(A_agg, F_agg, B_agg, C_agg, D_agg, E_agg, G_agg, H_agg, I_agg)
  ),

  tar_target(
    all_rep,
    dplyr::bind_rows(A_rep, F_rep, B_rep, C_rep, D_rep, E_rep, G_rep, H_rep, I_rep)
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

  # ---- Long-format scenario × model settings table ----
  # One row per (scenario_id, model_label) actually executed by the pipeline,
  # with one column per DGP and fit configuration field. Used to cross-
  # reference aggregated_results.csv against the exact config that produced
  # each row.
  tar_target(
    scenario_settings_csv,
    {
      dir.create("output", showWarnings = FALSE)
      tbl <- build_scenario_settings_table()
      readr::write_csv(tbl, "output/scenario_settings.csv")
      "output/scenario_settings.csv"
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
