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

# Individual-level data fits are markedly slower than summary fits (the
# per-subject likelihood plus per-study rho/effect parameters), so the
# scenarios that use them (A11, A12, I7) run at a reduced replication count.
# Categories A and I are split into a "main" map (N_REPS) and an "_indiv" map
# (N_REPS_INDIV), then recombined under the original A_rep / I_rep name so all
# downstream targets are unchanged.
N_REPS_INDIV <- 15L

# Figure-sweep categories (J-P) run at a deliberately small replication
# count while the paper figure design settles. Overridable via the
# N_REPS_FIG environment variable (used by the figure.yaml workflow's
# reps input); raise for the final run.
N_REPS_FIG <- as.integer(Sys.getenv("N_REPS_FIG", "5"))

# Category M feeds panel F, which plots an RMSE RATIO -- a second-moment
# statistic whose Monte Carlo error is ~1/sqrt(2n): about 32% at 5 reps,
# against ~8% for the posterior-SD and interval-width panels, which are means.
# At 5 reps that panel is unreadable while the others are fine.
#
# M is also the cheapest category to run (~2.2 min of fitting per rep vs ~7.2
# for Q), so it can take roughly 3x the reps and still finish inside Q's
# shadow. 3x is chosen as the largest multiplier that keeps M OFF the critical
# path at any global rep count -- so raising N_REPS_FIG for the camera-ready
# run does not make M the binding constraint or push it into the 300-minute
# job timeout ahead of Q.
N_REPS_FIG_M <- 3L * N_REPS_FIG

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
tar_source("R/illustration.R")

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
  # Summary-data scenarios at full reps; individual-data scenarios (A11, A12)
  # at reduced reps. Recombined into A_rep so downstream targets are unchanged.
  tarchetypes::tar_map_rep(
    name    = A_rep_main,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(function(s) !scenario_is_individual(s), scenario_ids("A"))),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tarchetypes::tar_map_rep(
    name    = A_rep_indiv,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(scenario_is_individual, scenario_ids("A"))),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_INDIV,
    reps    = 1
  ),
  tar_target(A_rep, dplyr::bind_rows(A_rep_main, A_rep_indiv)),
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
  # Excludes scenarios with correlated_effects = TRUE (incompatible with robust).
  # Split main/indiv so individual-data scenarios run at the reduced rep count,
  # matching their base (normal) arm; recombined into A_rep_robust.
  tarchetypes::tar_map_rep(
    name    = A_rep_robust_main,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(
      function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects) &&
        !scenario_is_individual(s),
      scenario_ids("A")
    )),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tarchetypes::tar_map_rep(
    name    = A_rep_robust_indiv,
    command = run_robust_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(
      function(s) !isTRUE(SCENARIO_CONFIGS[[s]]$fit$correlated_effects) &&
        scenario_is_individual(s),
      scenario_ids("A")
    )),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_INDIV,
    reps    = 1
  ),
  tar_target(A_rep_robust, dplyr::bind_rows(A_rep_robust_main, A_rep_robust_indiv)),
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
  # Summary-data I scenarios at full reps; individual-data I7 at reduced reps.
  # Recombined into I_rep so downstream targets are unchanged.
  tarchetypes::tar_map_rep(
    name    = I_rep_main,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(function(s) !scenario_is_individual(s), scenario_ids("I"))),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS,
    reps    = 1
  ),
  tarchetypes::tar_map_rep(
    name    = I_rep_indiv,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(Filter(scenario_is_individual, scenario_ids("I"))),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_INDIV,
    reps    = 1
  ),
  tar_target(I_rep, dplyr::bind_rows(I_rep_main, I_rep_indiv)),
  tar_target(I_agg, aggregate_scenario(I_rep)),

  # ---- Category J: figure sweep (panel B; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = J_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("J")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(J_agg, aggregate_scenario(J_rep)),

  # ---- Category K: figure sweep (panel D; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = K_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("K")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(K_agg, aggregate_scenario(K_rep)),

  # ---- Category L: figure sweep (panel E; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = L_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("L")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(L_agg, aggregate_scenario(L_rep)),

  # ---- Category M: figure sweep (panels C and F; extra reps) ----
  tarchetypes::tar_map_rep(
    name    = M_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("M")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG_M,
    reps    = 1
  ),
  tar_target(M_agg, aggregate_scenario(M_rep)),

  # ---- Category N: figure sweep (panel F; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = N_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("N")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(N_agg, aggregate_scenario(N_rep)),

  # ---- Category O: figure sweep (panel G; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = O_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("O")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(O_agg, aggregate_scenario(O_rep)),

  # ---- Category P: figure sweep (panel H; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = P_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("P")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(P_agg, aggregate_scenario(P_rep)),

  # ---- Category Q: trend-plane grid (panels K/L; reduced reps) ----
  tarchetypes::tar_map_rep(
    name    = Q_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("Q")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(Q_agg, aggregate_scenario(Q_rep)),

  # ---- Category R: DiD-count sweep at the pure null (panel J) ----
  tarchetypes::tar_map_rep(
    name    = R_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("R")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(R_agg, aggregate_scenario(R_rep)),

  # ---- Category U: per-design exchange rate vs rho (panel E) ----
  tarchetypes::tar_map_rep(
    name    = U_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("U")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(U_agg, aggregate_scenario(U_rep)),

  # ---- Category T: M sweep at a PP-heavy composition (panel F) ----
  tarchetypes::tar_map_rep(
    name    = T_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("T")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(T_agg, aggregate_scenario(T_rep)),

  # ---- Paired metadid/naive RMSE contrast for panel F ----
  # Computed from rep-level data, not from the aggregates: both arms share a
  # simulated dataset within a replicate, and resampling replicates keeps that
  # pairing (see paired_rmse_ratio()). M only -- referencing T_rep here would
  # silently pull the whole T category into any build of this target.
  tar_target(paired_rmse_M,  paired_ratio(M_rep, "bias", "rmse")),
  # Panel C is a ratio of mean interval widths -- a first-moment statistic, so
  # its interval is far tighter than panel F's RMSE ratio at the same reps.
  tar_target(paired_width_R, paired_ratio(R_rep, "ci_width", "mean")),
  tar_target(
    paired_contrasts,
    dplyr::bind_rows(paired_rmse_M, paired_width_R)
  ),

  # ---- Category S: composition sweep, high effect heterogeneity ----
  tarchetypes::tar_map_rep(
    name    = S_rep,
    command = run_one_rep(scenario_id, config, targets::tar_seed_get(), metadid_src),
    values  = scenario_values(scenario_ids("S")),
    names   = tidyselect::any_of("scenario_id"),
    batches = N_REPS_FIG,
    reps    = 1
  ),
  tar_target(S_agg, aggregate_scenario(S_rep)),

  # ---- One-dataset illustration (Figure 1, column 1) ----
  # Arguments NAMED deliberately: `pkg` exists only to create a dependency
  # edge on the metadid source, and a positional call put it in whichever
  # parameter happened to be second -- which broke when run_illustration()
  # gained n_studies/sigma_effect.
  tar_target(
    illustration_draws,
    run_illustration(seed = 495L, pkg = metadid_src)
  ),
  tar_target(
    illustration_csv,
    {
      dir.create("output", showWarnings = FALSE)
      readr::write_csv(illustration_draws, "output/illustration_draws.csv")
      "output/illustration_draws.csv"
    },
    format = "file"
  ),

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
  tar_target(diag_plot_J, plot_diagnostics(J_rep,  J_agg,  scenario_lookup_tbl, "J")),
  tar_target(diag_plot_K, plot_diagnostics(K_rep,  K_agg,  scenario_lookup_tbl, "K")),
  tar_target(diag_plot_L, plot_diagnostics(L_rep,  L_agg,  scenario_lookup_tbl, "L")),
  tar_target(diag_plot_M, plot_diagnostics(M_rep,  M_agg,  scenario_lookup_tbl, "M")),
  tar_target(diag_plot_N, plot_diagnostics(N_rep,  N_agg,  scenario_lookup_tbl, "N")),
  tar_target(diag_plot_O, plot_diagnostics(O_rep,  O_agg,  scenario_lookup_tbl, "O")),
  tar_target(diag_plot_P, plot_diagnostics(P_rep,  P_agg,  scenario_lookup_tbl, "P")),

  # ---- Combined results ----
  tar_target(
    all_agg,
    # Q/R/S are figure-only categories: they are generated programmatically
    # so sims.yaml's source grep cannot see them, which means no store slice
    # exists and the report job would re-simulate all 52 of them serially on
    # one runner. They are built by figure.yaml instead.
    dplyr::bind_rows(A_agg, F_agg, B_agg, C_agg, D_agg, E_agg, G_agg, H_agg,
                     I_agg, J_agg, K_agg, L_agg, M_agg, N_agg, O_agg, P_agg)
  ),

  tar_target(
    all_rep,
    dplyr::bind_rows(A_rep, F_rep, B_rep, C_rep, D_rep, E_rep, G_rep, H_rep,
                     I_rep, J_rep, K_rep, L_rep, M_rep, N_rep, O_rep, P_rep)
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

  # ---- Figure-panel data export ----
  # One tidy CSV with the treatment_effect_mean rows that feed the paper's
  # six-panel simulation figure (panel A from the existing A/F calibration
  # scenarios; panels B-F from the J-N sweeps).
  tar_target(
    figure_data_csv,
    {
      dir.create("output", showWarnings = FALSE)
      # Tagged by CATEGORY, not by panel position. Panel letters have already
      # been reshuffled once as panels were added and cut; category letters are
      # stable identity, so the figure code owns layout and this CSV does not
      # go stale when a panel moves. Q/R/S are figure-only (see all_agg).
      panels <- dplyr::bind_rows(
        dplyr::mutate(dplyr::bind_rows(A_agg, F_agg), category = "AF"),
        dplyr::mutate(J_agg, category = "J"),
        dplyr::mutate(M_agg, category = "M"),
        dplyr::mutate(K_agg, category = "K"),
        dplyr::mutate(L_agg, category = "L"),
        dplyr::mutate(N_agg, category = "N"),
        dplyr::mutate(O_agg, category = "O"),
        dplyr::mutate(P_agg, category = "P")
      ) |>
        dplyr::filter(parameter == "treatment_effect_mean")
      readr::write_csv(panels, "output/figure_panels.csv")
      "output/figure_panels.csv"
    },
    format = "file"
  ),

  # ---- Validation report ----
  tarchetypes::tar_quarto(
    report,
    path = "reports/validation-report.qmd"
  )
)
