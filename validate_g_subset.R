# Focused validation of the G-series ("Jensen") scenarios after correcting the
# normalised truth to the per-study percentage estimand E[theta_i/b_i].
#
# For each scenario it fits the (per-study normalised) metadid model over a few
# replications and reports the mean bias of treatment_effect_mean against
#   * the CORRECTED truth  E[theta_i/b_i]   (mean of per-study % effects), and
#   * the OLD truth        E[theta]/E[b]     (= true_effect / baseline_mean).
# The difference between the two biases is the deterministic between-study
# Jensen term, true_effect/bm * CV_b^2. Expectation: bias_new ~ 0 at every
# scenario; bias_old ~ 0 at G6 (baseline_sd = 0) and ~ -0.016 at G8
# (baseline_sd = 0.10) -- i.e. the historical "G8 bias" was scoring against the
# wrong target. See ben18785/metadid#39.
#
# Requires the corrected (per-study) metadid installed and the sibling ../metadid
# on the revert/correct-estimand branch. Run from the metadid-sims checkout:
#   Rscript validate_g_subset.R
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(metadid)
})
# Source only the files needed to simulate + fit (avoid R/plots.R, which needs
# kableExtra / a full reporting toolchain).
for (f in c("R/scenarios.R", "R/simulate.R", "R/fit.R")) source(f)
`%or%` <- function(a, b) if (is.null(a)) b else a

run_scenario <- function(sc, nreps) {
  cfg       <- SCENARIO_CONFIGS[[sc]]
  bm        <- cfg$dgp$baseline_mean
  old_truth <- cfg$dgp$true_effect / bm  # E[theta]/E[b] (previous, wrong target)
  rows <- map_dfr(seq_len(nreps), function(r) {
    sim  <- simulate_scenario(sc, rep_seed = 7000L + r, cfg)
    post <- fit_scenario(sim, cfg$fit)
    te   <- post[post$parameter == "treatment_effect_mean", ]
    tibble(
      new_truth = sim$true_params$treatment_effect_mean_normalised,  # corrected
      post_mean = te$mean, lo = te$lo, hi = te$hi
    )
  })
  tibble(
    scenario    = sc,
    baseline_sd = cfg$dgp$baseline_sd %or% 0.05,
    n_did       = cfg$dgp$n_did %or% 20L,
    reps        = nreps,
    new_truth   = round(mean(rows$new_truth), 4),
    old_truth   = round(old_truth, 4),
    bias_new    = round(mean(rows$post_mean - rows$new_truth), 4),
    bias_old    = round(mean(rows$post_mean - old_truth), 4),
    cov_new     = round(mean(rows$lo < rows$new_truth & rows$hi > rows$new_truth), 2),
    cov_old     = round(mean(rows$lo < old_truth      & rows$hi > old_truth),      2)
  )
}

reps <- list(G6 = 6L, G7 = 6L, G8 = 10L)
res  <- bind_rows(lapply(names(reps), function(s) run_scenario(s, reps[[s]])))
cat("\n===== G-series subset: per-study model vs corrected (new) and old truths =====\n")
print(as.data.frame(res), row.names = FALSE)
cat("\nDONE\n")
