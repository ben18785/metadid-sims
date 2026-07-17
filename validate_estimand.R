# Focused estimand validation for the baseline-weighted reporting fix.
#
# Reproduces the metadid-sims "Jensen" G-scenarios (between-study baseline
# variation) on a handful of replications. The default per-study normalised
# fit now REPORTS the baseline-weighted mean of per-study fractional effects,
# i.e. the ratio of population means E[theta]/E[b] (summary()'s
# treatment_effect_mean), which is the unbiased target estimand. The raw Stan
# parameter treatment_effect_mean still estimates the equal-weighted mean
# E[theta_i/b_i], which carries the Jensen bias. Truth = -0.15/0.45 = -0.3333.
#
# Expectation (per replication, averaged over reps):
#   * G6 (baseline_sd = 0):    weighted ~ raw ~ unbiased (no Jensen term).
#   * G8 (baseline_sd = 0.10): weighted ~unbiased; raw biased low by ~ -0.016.
#   * Convergence is the ordinary per-study normalised fit (robust): Rhat ~1.0
#     on every replication (the reweighting is pure post-processing).
#
# Run from the metadid-sims checkout with sibling ../metadid on the
# fix-weighted-estimand branch:  Rscript validate_estimand.R
suppressMessages({library(devtools); load_all("../metadid", quiet = TRUE)})

TRUTH <- -0.15 / 0.45
NREPS <- 10
SCEN  <- c(G6 = 0.0, G8 = 0.10)

run_cell <- function(bsd) {
  wbias <- rbias <- cover <- rhat <- numeric(NREPS)
  for (r in seq_len(NREPS)) {
    sim <- simulate_meta_did(
      n_studies = 100, true_effect = -0.15, sigma_effect = 0.03,
      true_trend = -0.04, sigma_trend = 0.02, baseline_mean = 0.45,
      baseline_sd = bsd, n_control = 100L, n_treatment = 100L, seed = 1000L + r
    )
    fit <- meta_did(summary_data = as_summary_did(sim), chains = 2L,
                    iter_warmup = 400L, iter_sampling = 400L, seed = 99L, refresh = 0)
    ss  <- summary(fit); te <- ss[ss$parameter == "treatment_effect_mean", ]
    raw <- mean(fit$fit$draws("treatment_effect_mean", format = "matrix"))
    wbias[r] <- te$mean - TRUTH
    rbias[r] <- raw - TRUTH
    cover[r] <- te$lo < TRUTH && te$hi > TRUTH
    rhat[r]  <- max(fit$fit$summary("treatment_effect_mean")$rhat)
  }
  sprintf("weighted_bias=%+.4f (cover %.2f) | raw_bias=%+.4f | max_rhat=%.3f",
          mean(wbias), mean(cover), mean(rbias), max(rhat))
}

cat(sprintf("TRUTH E[theta]/E[b] = %.4f ; %d reps/cell\n\n", TRUTH, NREPS))
for (nm in names(SCEN)) {
  cat(sprintf("%s (baseline_sd=%.2f): %s\n", nm, SCEN[[nm]], run_cell(SCEN[[nm]])))
}
cat("\nDONE\n")
