# ===========================================================================
# Diagnostic: where does the post-only residual bias come from?
# ===========================================================================
#
# The validation report carries one open issue: a bias of roughly +0.02 on a
# true effect of -0.337 (about 7%, toward zero) whenever post-only studies
# dominate the evidence base. It is absent from DiD-only scenarios and from
# category A, and the delta-method normalisation patch did not move it.
#
# The leading suspicion is the normalisation denominator: post-only studies
# divide by mean_post_control (alpha + beta) while DiD studies divide by
# mean_pre_control (alpha). The model is meant to undo that mismatch with
# theta_tilde = phi * (1 + beta_tilde), so a mis-applied factor there would
# show up exactly as a bias in the mean.
#
# Two facts already narrow this down, from the published aggregate CSV:
#
#   * time_trend_mean is NOT biased in the affected arms (-0.005 to +0.014
#     against a true -0.089, coverage 0.75-1.00). Producing +0.02 through a
#     mis-estimated trend would need the trend wrong by about -0.054, so
#     "the trend is wrong and the conversion amplifies it" is ruled out.
#
#   * The arithmetic brackets the observed value. Writing th = theta_tilde and
#     b = beta_tilde:
#         conversion omitted  -> th/(1+b) - th      = -0.033  (wrong sign)
#         conversion correct  ->                       0
#         conversion doubled  -> th*(1+b) - th      = +0.030  (observed sign)
#
# This script tests the prediction that distinguishes a mis-applied conversion
# factor from everything else: the bias must be PROPORTIONAL TO THE TRUE TREND,
# and must FLIP SIGN when the trend does. Nothing in the gamma machinery or in
# precision weighting predicts that.
#
# Note on magnitudes: the predictions printed below are the full-strength
# effect on a post-only study. The pooled estimate blends those studies with
# unaffected DiD studies, so the realised bias is a fraction of the printed
# value (roughly the post-only share of the evidence). The LINEARITY IN THE
# TREND is the diagnostic, not the absolute size.
#
# Composition is 10 DiD + 20 post-only throughout, using the standard metadid
# DGP, where gamma is identically zero for every study. That isolates the
# normalisation conversion from the baseline-imbalance machinery: the reported
# issue is present even in arms whose post-only studies have gamma == 0.
#
# Usage:  Rscript diagnostics/post-only-bias.R
#   REPS        replications per arm (default 30)
#   DIAG_CORES  concurrent fits (default 2; each fit itself uses 2 chains)
#   OUT_DIR     where to write outputs (default "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

REPS       <- as.integer(Sys.getenv("REPS", "30"))
DIAG_CORES <- as.integer(Sys.getenv("DIAG_CORES", "2"))
OUT_DIR    <- Sys.getenv("OUT_DIR", "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_DID <- 10L
N_RCT <- 20L

# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------
# T* sweep the true trend with the default fit. F* vary the fit at the default
# trend. C* are harness controls that should show no bias either way.

ARMS <- list(

  T_trend_0 = scenario(
    "Trend 0: 10 DiD + 20 post-only, no time trend at all",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = 0, sigma_trend = 0)
  ),

  T_trend_m04 = scenario(
    "Trend -0.04 (pipeline default): the arm that should reproduce the bias",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04)
  ),

  T_trend_m08 = scenario(
    "Trend -0.08: double the default, so a conversion fault should double",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.08)
  ),

  T_trend_p04 = scenario(
    "Trend +0.04: the sign-flip test, the decisive arm",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = +0.04)
  ),

  # The tau_beta axis. A post-only study's own data say nothing about beta, so
  # its conversion factor is inherited: SD(beta_tilde_i) ~ sqrt(tau_beta^2 +
  # SE(mu_beta)^2). Shrinking tau_beta shrinks both terms at once -- the spread
  # directly, and mu_beta because the DiD studies' precision goes as
  # n_did / (s_beta^2 + tau_beta^2). At tau_beta = 0 every study's trend IS the
  # population trend, so borrowing stops being an approximation and the factor
  # is effectively known.
  #
  # This separates two faults that the trend sweep alone cannot, because both
  # scale with beta: a STRUCTURAL one (an extra factor applied where it should
  # not be) depends on beta's value and survives tau_beta = 0, whereas an
  # UNCERTAINTY-DRIVEN one (posterior geometry tilting a diffuse factor) does
  # not. T_trend_m04 is the midpoint of this axis (sigma_trend = 0.02).

  H_tau_b_0 = scenario(
    "No trend heterogeneity: the conversion factor is effectively known",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04, sigma_trend = 0)
  ),

  H_tau_b_04 = scenario(
    "Double trend heterogeneity: a more diffuse conversion factor",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04, sigma_trend = 0.04)
  ),

  H_more_did = scenario(
    "30 DiD + 20 post-only: sharpens mu_beta without touching tau_beta",
    dgp = list(n_did = 30L, n_rct = N_RCT, true_trend = -0.04)
  ),

  F_raw = scenario(
    "Unnormalised fit at the default trend: no denominator, so no conversion",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04),
    fit = list(normalise = FALSE)
  ),

  # fn MUST be switched here: fit_scenario() only forwards time_trend when
  # fn == "meta_did_general" (R/fit.R), so leaving the default "meta_did" makes
  # the override a silent no-op. The 3-rep smoke run caught this -- the arm came
  # back byte-identical to T_trend_m04.
  F_fixed_zero = scenario(
    "time_trend = fixed_zero: bypasses the reparameterisation entirely",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04),
    fit = list(fn = "meta_did_general", time_trend = "fixed_zero")
  ),

  # Is there a bias floor common to every arm, independent of design? Removing
  # baseline variation removes all spread in the normalisation denominators, so
  # both the Jensen correction to the estimand and the precision-weighting
  # effect (the mechanism that moved Fit 6) vanish. Whatever survives here is
  # neither of those.
  Z_bsd0 = scenario(
    "baseline_sd = 0: no denominator spread, so no weighting floor",
    dgp = list(n_did = N_DID, n_rct = N_RCT, true_trend = -0.04, baseline_sd = 0)
  ),

  # The arms above all use the STANDARD metadid DGP. The published +0.02 lives
  # in X1-X9 and X13-X18, which use the bespoke simulate_imbalance_grid. These
  # three mirror that generator so we can tell whether the phenomenon is even in
  # scope here. kappa_true = 0 makes gamma identically zero for the randomised
  # studies, and the did_gamma_* defaults are zero too, so no study carries any
  # baseline imbalance -- which is the configuration the open issue says still
  # shows the bias.
  B_grid = scenario(
    "Bespoke imbalance-grid DGP, gamma == 0 everywhere: does +0.02 reproduce?",
    dgp = list(type = "bespoke", bespoke_fn = "simulate_imbalance_grid",
               n_did = N_DID, n_rct = N_RCT, n_randomised_did = 0L,
               n_randomised_rct = N_RCT, kappa_true = 0, true_trend = -0.04)
  ),

  B_grid_p04 = scenario(
    "As B_grid with trend +0.04: the sign-flip test on the DGP that shows it",
    dgp = list(type = "bespoke", bespoke_fn = "simulate_imbalance_grid",
               n_did = N_DID, n_rct = N_RCT, n_randomised_did = 0L,
               n_randomised_rct = N_RCT, kappa_true = 0, true_trend = +0.04)
  ),

  B_grid_bsd0 = scenario(
    "As B_grid with baseline_sd = 0: bespoke DGP without the weighting floor",
    dgp = list(type = "bespoke", bespoke_fn = "simulate_imbalance_grid",
               n_did = N_DID, n_rct = N_RCT, n_randomised_did = 0L,
               n_randomised_rct = N_RCT, kappa_true = 0, true_trend = -0.04,
               baseline_sd = 0)
  ),

  C_did_only = scenario(
    "Control: 30 DiD, no post-only studies — expect no bias",
    dgp = list(n_did = 30L, n_rct = 0L, true_trend = -0.04)
  )
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

run_arm <- function(arm_id, config, reps) {
  message(sprintf("[%s] %s", arm_id, config$description))
  out <- mclapply(seq_len(reps), function(r) {
    res <- try(run_one_rep(arm_id, config, rep_seed = 20261029L + r), silent = TRUE)
    if (inherits(res, "try-error")) {
      message(sprintf("  rep %d failed: %s", r, conditionMessage(attr(res, "condition"))))
      return(NULL)
    }
    res$rep <- r
    res
  }, mc.cores = max(1L, DIAG_CORES))
  ok <- out[!vapply(out, is.null, logical(1))]
  if (length(ok) == 0) {
    message(sprintf("  [%s] every replication failed", arm_id))
    return(NULL)
  }
  if (length(ok) < reps)
    message(sprintf("  [%s] %d of %d replications usable", arm_id, length(ok), reps))
  bind_rows(ok)
}

reps_all <- bind_rows(lapply(names(ARMS), function(id) run_arm(id, ARMS[[id]], REPS)))
stopifnot(nrow(reps_all) > 0)
agg <- aggregate_scenario(reps_all)

# ---------------------------------------------------------------------------
# Predictions
# ---------------------------------------------------------------------------
# For a normalised fit, th = true_value of treatment_effect_mean and
# b = true_trend / baseline_mean. The unnormalised arm has no conversion, and
# the fixed_zero arm declines to apply one, so each gets a single prediction.

arm_meta <- bind_rows(lapply(names(ARMS), function(id) {
  d <- ARMS[[id]]$dgp
  tibble(scenario_id = id, description = ARMS[[id]]$description,
         normalise = ARMS[[id]]$fit$normalise,
         time_trend = ARMS[[id]]$fit$time_trend,
         n_did = d$n_did, n_rct = d$n_rct,
         true_trend = d$true_trend, sigma_trend = d$sigma_trend,
         b_norm = d$true_trend / d$baseline_mean)
}))

effect <- agg |>
  filter(parameter == "treatment_effect_mean") |>
  select(scenario_id, n_reps, true_value, empirical_coverage,
         mean_bias, se_bias, rmse, pct_rhat_bad, mean_n_divergent) |>
  left_join(arm_meta, by = "scenario_id") |>
  mutate(
    # Full-strength predictions for a post-only study; the pooled bias is a
    # fraction of these (the post-only share of the evidence).
    pred_doubled = if_else(normalise & time_trend != "fixed_zero",
                           true_value * b_norm, NA_real_),
    pred_omitted = if_else(normalise,
                           -true_value * b_norm / (1 + b_norm), NA_real_),
    bias_over_b  = if_else(b_norm != 0, mean_bias / b_norm, NA_real_)
  )

trend_rows <- agg |>
  filter(parameter == "time_trend_mean") |>
  select(scenario_id, trend_true = true_value,
         trend_bias = mean_bias, trend_cov = empirical_coverage)

tau_rows <- agg |>
  filter(parameter == "treatment_effect_sd") |>
  select(scenario_id, tau_true = true_value, tau_bias = mean_bias)

summary_tbl <- effect |>
  left_join(trend_rows, by = "scenario_id") |>
  left_join(tau_rows, by = "scenario_id") |>
  arrange(match(scenario_id, names(ARMS)))

readr::write_csv(summary_tbl, file.path(OUT_DIR, "post-only-bias-diagnostic.csv"))
readr::write_csv(reps_all,    file.path(OUT_DIR, "post-only-bias-reps.csv"))

# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

fmt <- function(x, d = 4) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))

md <- c(
  "## Post-only bias diagnostic",
  "",
  sprintf("%d replications per arm, %d DiD + %d post-only studies unless noted.",
          REPS, N_DID, N_RCT),
  "",
  "`bias/b` is the observed bias divided by the true normalised trend. If a",
  "mis-applied conversion factor is the cause, this column is roughly CONSTANT",
  "across the trend sweep (and the +0.04 arm's bias is negative).",
  "",
  "| Arm | DiD/post | trend | sigma_b | reps | bias | MC se | bias/b | pred if doubled | pred if omitted | coverage | trend bias |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d/%d | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |",
          summary_tbl$scenario_id, summary_tbl$n_did, summary_tbl$n_rct,
          fmt(summary_tbl$true_trend, 2), fmt(summary_tbl$sigma_trend, 2),
          summary_tbl$n_reps,
          fmt(summary_tbl$mean_bias), fmt(summary_tbl$se_bias),
          fmt(summary_tbl$bias_over_b, 3),
          fmt(summary_tbl$pred_doubled), fmt(summary_tbl$pred_omitted),
          fmt(summary_tbl$empirical_coverage, 2),
          fmt(summary_tbl$trend_bias)),
  "",
  "### How to read it",
  "",
  "* **`T_trend_p04` bias is negative** — the conversion factor is the cause.",
  "  Compare `bias/b` across the four `T_` arms: near-constant means the fault",
  "  scales with the trend exactly as a mis-applied `(1 + beta_tilde)` would.",
  "* **`T_trend_0` bias is near zero while the others are not** — same verdict,",
  "  since with no trend there is nothing for the conversion to get wrong.",
  "* **Bias stays positive across all four `T_` arms** — the conversion is",
  "  exonerated; look to the gamma machinery and to precision weighting instead.",
  "* **`F_raw` is clean** — consistent with normalisation being the locus",
  "  (compare X24-X27, whose raw arms sit at 0.008 against 0.024-0.043).",
  "* **`F_fixed_zero` lands near `pred if omitted`** — the conversion really is",
  "  needed and the normalised data does measure `theta/(alpha + beta)`. If it",
  "  comes out near zero instead, the denominator is not what we think it is and",
  "  the fault is upstream, in the normalisation path rather than the model.",
  "* **`C_did_only` must be near zero.** If it is not, the harness is at fault",
  "  and nothing else here can be trusted.",
  "",
  "### The tau_beta axis: structural or uncertainty-driven?",
  "",
  "Read `H_tau_b_0` (sigma_b = 0), `T_trend_m04` (0.02) and `H_tau_b_04` (0.04)",
  "together, all at the same true trend. `H_more_did` sharpens mu_beta instead,",
  "leaving tau_beta alone, and should move the same way as `H_tau_b_0` if the",
  "fault is uncertainty-driven.",
  "",
  "* **Bias survives sigma_b = 0** -> STRUCTURAL. The conversion factor's value",
  "  is wrong, not merely uncertain, and the fault is in the mean structure.",
  "  Combined with a bias that is linear in the trend, that points at an extra",
  "  `(1 + beta_tilde)` applied somewhere it should not be.",
  "* **Bias shrinks toward zero as sigma_b does** -> UNCERTAINTY-DRIVEN. The",
  "  posterior geometry is tilting a diffuse borrowed factor, and the remedy is",
  "  a different one (parameterisation, not arithmetic).",
  "* **Bias is flat in sigma_b AND flat in the trend** -> neither; the conversion",
  "  is exonerated and precision weighting becomes the leading suspect. That",
  "  mechanism is already confirmed to exist (it is what moved Fit 6 and the",
  "  X1-X6 arms after the delta-method patch), but at baseline_sd = 0.05 it",
  "  accounts for only about +0.004 of the roughly +0.02 observed.",
  "",
  "### Is there a floor, and is the phenomenon even in scope?",
  "",
  "`C_did_only` and `Z_bsd0` test for a bias common to every arm. If all arms",
  "sit at the same small positive value, that is a FLOOR rather than a post-only",
  "effect, and `Z_bsd0` says what it is: with no spread in the denominators,",
  "both the Jensen correction and the weighting effect are gone, so anything",
  "left over is a third thing.",
  "",
  "`B_grid*` mirror simulate_imbalance_grid, the generator the published +0.02",
  "arms actually use, with gamma identically zero throughout. They answer a",
  "prior question: does the standard DGP reproduce the phenomenon at all?",
  "",
  "* **B_grid near +0.02 while the standard-DGP arms sit near the floor** -> the",
  "  bias is specific to that generator, and the whole trend/tau_beta sweep above",
  "  was run on data that never exhibited it. Rebuild the sweep on the bespoke",
  "  DGP before concluding anything from it.",
  "* **B_grid at the same floor as everything else** -> the +0.02 in the report",
  "  comes from something those scenarios add beyond the generator (non-zero",
  "  gamma, randomisation labels, kappa_true > 0, cluster settings), and the next",
  "  step is to switch those on one at a time.",
  "* **B_grid_p04 flips sign while B_grid does not** -> the conversion is",
  "  implicated after all, but only under this generator.",
  "",
  sprintf("Convergence: worst %%Rhat>1.05 = %.2f, worst mean divergences = %.1f.",
          max(summary_tbl$pct_rhat_bad, na.rm = TRUE),
          max(summary_tbl$mean_n_divergent, na.rm = TRUE)),
  "",
  "Monte Carlo error matters: at 30 replications the MC se column is the",
  "yardstick, and differences smaller than about twice it are noise."
)

writeLines(md, file.path(OUT_DIR, "post-only-bias-diagnostic.md"))
cat(paste(md, collapse = "\n"), "\n")
