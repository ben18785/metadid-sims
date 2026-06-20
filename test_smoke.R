# test_smoke.R — fast smoke test for the simulation pipeline
#
# Runs one replication of a representative scenario from each category and
# function variant (run_one_rep, run_robust_rep, run_naive_rep, run_g_rep)
# using minimal MCMC settings. Takes ~5-10 minutes total rather than hours.
#
# Usage:
#   source("test_smoke.R")

library(targets)
library(metadid)

tar_source("R/scenarios.R")
tar_source("R/simulate.R")
tar_source("R/fit.R")
tar_source("R/fit_g.R")
tar_source("R/fit_robust.R")
tar_source("R/fit_naive.R")
tar_source("R/assess.R")

# Override MCMC settings for speed — do this AFTER sourcing fit.R
MCMC_OPTS <<- list(
  chains          = 2L,
  parallel_chains = 2L,
  iter_warmup     = 200L,
  iter_sampling   = 200L,
  refresh         = 0
)

# ---------------------------------------------------------------------------
# Helper: run one test case and report pass/fail
# ---------------------------------------------------------------------------

run_test <- function(label, expr) {
  cat(sprintf("  %-55s", label))
  result <- tryCatch(
    { force(expr); "PASS" },
    error   = function(e) paste("FAIL:", conditionMessage(e)),
    warning = function(w) "PASS (with warning)"
  )
  cat(result, "\n")
  invisible(result)
}

cat("\n=== metadid-sims smoke test ===\n\n")

# ---------------------------------------------------------------------------
# One scenario per category via run_one_rep
# ---------------------------------------------------------------------------

cat("run_one_rep (standard scenarios):\n")
for (sid in c("A1", "B1", "C1", "D1", "E1", "F1", "H1", "H2", "H3", "H4", "I1", "I8")) {
  cfg <- SCENARIO_CONFIGS[[sid]]
  run_test(sid, run_one_rep(sid, cfg, rep_seed = 1L))
}

# ---------------------------------------------------------------------------
# Robust comparator
# ---------------------------------------------------------------------------

cat("\nrun_robust_rep:\n")
for (sid in c("A1", "F1")) {
  cfg <- SCENARIO_CONFIGS[[sid]]
  run_test(sid, run_robust_rep(sid, cfg, rep_seed = 1L))
}

# ---------------------------------------------------------------------------
# Naive comparator
# ---------------------------------------------------------------------------

cat("\nrun_naive_rep:\n")
for (sid in c("A1", "A9")) {  # A9 is unnormalised
  cfg <- SCENARIO_CONFIGS[[sid]]
  run_test(sid, run_naive_rep(sid, cfg, rep_seed = 1L))
}

# ---------------------------------------------------------------------------
# G scenarios (prior sensitivity + Jensen's inequality)
# ---------------------------------------------------------------------------

cat("\nrun_g_rep:\n")
for (sid in c("G1", "G3", "G6", "G9")) {
  cfg <- SCENARIO_CONFIGS[[sid]]
  run_test(sid, run_g_rep(sid, cfg, rep_seed = 1L))
}

cat("\nDone.\n")
