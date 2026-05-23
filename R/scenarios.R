# Scenario definitions for the metadid simulation study.
#
# Each scenario is a named list stored in SCENARIO_CONFIGS. The structure:
#
#   dgp:  Data-generating process parameters
#     type          "metadid" (use simulate_meta_did) or "bespoke"
#     bespoke_fn    Name of bespoke simulation function (when type = "bespoke")
#     n_did, n_rct, n_pp  Number of studies per design
#     true_effect, sigma_effect, true_trend, sigma_trend
#     baseline_mean, baseline_sd, within_sd, rho
#     rho_effect_trend
#     n_control, n_treatment
#     covariates, beta_cov    (optional)
#     ... any bespoke-specific params
#
#   fit:  Model fitting options
#     fn              "meta_did" or "meta_did_general"
#     normalise              logical — TRUE expresses effects as fractions of
#                            the treatment-pre baseline, FALSE pools on the
#                            absolute (user-units) scale
#     baseline_latent_arm    "treatment" (default) or "control" — which arm's
#                            baseline is the per-study latent in modelled mode
#     robust_heterogeneity, design_effects, correlated_effects
#     hierarchical_rho
#     time_trend, baseline_imbalance, pp_likelihood (for meta_did_general)
#     covariates      (formula or NULL)
#     provide_rho     Whether to include rho in summary data (default TRUE)
#     data_format     "summary" or "individual"
#
#   true:  Named list of true parameter values on the RAW scale.
#          The pipeline computes normalised values from these + baseline_mean.
#
#   compare:  (optional) For comparative studies, a list of alternative fit configs.
#             Each element is a named list like `fit` above, with a `label` field.

# ---------------------------------------------------------------------------
# Default DGP and fit values
# ---------------------------------------------------------------------------

default_dgp <- list(
  type             = "metadid",
  n_did            = 20L,
  n_rct            = 0L,
  n_pp             = 0L,
  true_effect      = -0.15,
  sigma_effect     = 0.03,
  true_trend       = -0.04,
  sigma_trend      = 0.02,
  baseline_mean    = 0.45,
  baseline_sd      = 0.05,
  within_sd        = 0.12,
  rho              = 0.5,
  rho_effect_trend = 0,
  n_control        = 100L,
  n_treatment      = 100L,
  covariates       = NULL,
  beta_cov         = NULL
)

default_fit <- list(
  fn                    = "meta_did",
  normalise             = TRUE,
  baseline_latent_arm   = "treatment",
  robust_heterogeneity  = FALSE,
  design_effects        = FALSE,
  correlated_effects    = FALSE,
  hierarchical_rho      = TRUE,
  time_trend            = "pooled",
  baseline_imbalance    = "estimated",
  pp_likelihood         = "differenced",
  covariates            = NULL,
  provide_rho           = TRUE,
  data_format           = "summary"
)

# Helper: merge user overrides into defaults
scenario <- function(description, dgp = list(), fit = list(),
                     true = list(), compare = NULL) {
  list(
    description = description,
    dgp         = modifyList(default_dgp, dgp),
    fit         = modifyList(default_fit, fit),
    true        = true,
    compare     = compare
  )
}

# ---------------------------------------------------------------------------
# Category A: Calibration studies
# ---------------------------------------------------------------------------

SCENARIO_CONFIGS <- list(

  A1 = scenario(
    "Baseline: DiD-only summary, normalised, 20 studies"
  ),

  A2 = scenario(
    "Large sample: DiD-only, 60 studies",
    dgp = list(n_did = 60L)
  ),

  A3 = scenario(
    "Small sample: DiD-only, 5 studies",
    dgp = list(n_did = 5L)
  ),

  A4 = scenario(
    "Very large sample: DiD-only, 100 studies",
    dgp = list(n_did = 100L)
  ),

  A5 = scenario(
    "Mixed designs: 10 DiD + 10 RCT + 10 PP, normalised",
    dgp = list(n_did = 10L, n_rct = 10L, n_pp = 10L)
  ),

  A6 = scenario(
    "PP-heavy with trend: 5 DiD + 25 PP, trend = -0.10",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.10, sigma_trend = 0.01)
  ),

  A7 = scenario(
    "Covariates: 30 DiD studies with dose covariate",
    dgp = list(
      n_did      = 30L,
      covariates = data.frame(dose = seq(1, 4, length.out = 30)),
      beta_cov   = -0.04
    ),
    fit = list(covariates = ~ dose)
  ),

  A8 = scenario(
    "Correlated effects: 30 DiD, rho_effect_trend = 0.95",
    dgp = list(n_did = 30L, rho_effect_trend = 0.95, baseline_sd = 0),
    fit = list(fn = "meta_did_general", correlated_effects = TRUE)
  ),

  A9 = scenario(
    "Unnormalised DiD-only, 20 studies, raw scale",
    fit = list(normalise = FALSE)
  ),

  A10 = scenario(
    "Unnormalised mixed: 10 DiD + 10 RCT + 10 PP, raw scale",
    dgp = list(n_did = 10L, n_rct = 10L, n_pp = 10L),
    fit = list(normalise = FALSE)
  ),

  A11 = scenario(
    "Individual-level DiD: fit from individual data, 20 studies",
    fit = list(data_format = "individual")
  ),

  A12 = scenario(
    "Individual vs summary consistency: same data, fit both ways, 20 studies",
    compare = list(
      list(label = "summary",    data_format = "summary"),
      list(label = "individual", data_format = "individual")
    )
  ),

  A13 = scenario(
    "Known vs hierarchical rho: same DiD data, compare rho modes",
    compare = list(
      list(label = "hierarchical_rho", hierarchical_rho = TRUE,  provide_rho = TRUE),
      list(label = "known_rho",        hierarchical_rho = FALSE, provide_rho = TRUE)
    )
  ),

  # ---------------------------------------------------------------------------
  # Category F: Large-N bias probes (separate targets pattern)
  # ---------------------------------------------------------------------------

  F1 = scenario(
    "Large mixed, normalised, full model: 35 DiD + 35 RCT + 30 PP",
    dgp = list(n_did = 35L, n_rct = 35L, n_pp = 30L)
  ),

  F2 = scenario(
    "Large mixed, normalised, naive: 70 DiD + 70 RCT + 60 PP (zero trend, equal baselines)",
    dgp = list(
      n_did = 70L, n_rct = 70L, n_pp = 60L,
      true_trend = 0, sigma_trend = 0, baseline_sd = 0
    ),
    fit = list(
      fn                 = "meta_did_general",
      time_trend         = "fixed_zero",
      baseline_imbalance = "fixed_zero"
    )
  ),

  F3 = scenario(
    "Large mixed, unnormalised: 70 DiD + 70 RCT + 60 PP",
    dgp = list(n_did = 70L, n_rct = 70L, n_pp = 60L),
    fit = list(normalise = FALSE)
  ),

  F4 = scenario(
    "Large DiD + RCT only: 50 DiD + 50 RCT",
    dgp = list(n_did = 50L, n_rct = 50L)
  ),

  F5 = scenario(
    "Large DiD + PP only: 100 DiD + 100 PP, nonzero trend",
    dgp = list(n_did = 100L, n_pp = 100L, true_trend = -0.10, sigma_trend = 0.01)
  ),

  F6 = scenario(
    "Large DiD-only with covariates: 200 studies",
    dgp = list(
      n_did      = 200L,
      covariates = data.frame(dose = seq(1, 4, length.out = 200)),
      beta_cov   = -0.04
    ),
    fit = list(covariates = ~ dose)
  ),

  F7 = scenario(
    "Large mixed, design effects: 35 DiD + 35 RCT + 30 PP with offsets",
    dgp = list(n_did = 35L, n_rct = 35L, n_pp = 30L),
    fit = list(design_effects = TRUE),
    true = list(delta_rct = 0.10, delta_pp = -0.08)
  ),

  F8 = scenario(
    "Large DiD-only, correlated effects: 200 studies, rho_effect_trend = 0.95",
    dgp = list(n_did = 200L, rho_effect_trend = 0.95, baseline_sd = 0),
    fit = list(fn = "meta_did_general", correlated_effects = TRUE)
  ),

  # ---------------------------------------------------------------------------
  # Category B: Comparative studies
  # ---------------------------------------------------------------------------

  B1 = scenario(
    "Naive vs full: PP-heavy, large trend (full should win)",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.10, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  B2 = scenario(
    "Naive vs full: assumptions satisfied (zero trend, equal baselines)",
    dgp = list(
      n_did = 10L, n_rct = 10L, n_pp = 10L,
      true_trend = 0, sigma_trend = 0, baseline_sd = 0
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  B3 = scenario(
    "Naive wins: DiD trend = 0, PP trend = -0.10 (pooled borrows wrong trend)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = 0, pp_trend = -0.10
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  B4 = scenario(
    "Correlated vs independent effects: rho_effect_trend = 0.95",
    dgp = list(n_did = 30L, rho_effect_trend = 0.95, baseline_sd = 0),
    compare = list(
      list(label = "correlated",  fn = "meta_did_general", correlated_effects = TRUE),
      list(label = "independent", fn = "meta_did")
    )
  ),

  B5 = scenario(
    "RCT baseline imbalance (only in this design): full (estimated) vs naive (fixed zero), 10 DiD + 10 RCT",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 10L,
      rct_gamma_mean = 0.05,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  B6 = scenario(
    "DiD imbalance > RCT imbalance: full vs naive, 15 DiD + 15 RCT",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_did_rct_imbalance",
      n_did          = 15L,
      n_rct          = 15L,
      did_gamma_mean = 0.08,
      did_gamma_sd   = 0.02,
      rct_gamma_mean = 0.01,
      rct_gamma_sd   = 0.01
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  # B7: paired comparison of the two modelled-mode parameterisations.
  # The two options of baseline_latent_arm encode the same statistical model
  # — they differ only in which arm's pre-baseline is the per-study latent
  # with the wide uniform prior, and which is derived via the hierarchical
  # baseline-imbalance parameter. In well-identified problems (decent sample
  # sizes, both arms observed) the two should give numerically equivalent
  # posteriors on every population-level parameter. This scenario fits the
  # same simulated data under both choices and lets the validation report
  # show the agreement (or any discrepancy under stress regimes).
  B7 = scenario(
    "Baseline-latent-arm parameterisation: treatment vs control, 30 DiD",
    dgp = list(n_did = 30L),
    compare = list(
      list(label = "treatment_latent", baseline_latent_arm = "treatment"),
      list(label = "control_latent",   baseline_latent_arm = "control")
    )
  ),

  # ---------------------------------------------------------------------------
  # Category C: Outlier and heavy-tailed studies
  # ---------------------------------------------------------------------------

  C1 = scenario(
    "Single outlier: one study shifted by 5*sigma_effect",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 1L,
      outlier_shift = 5   # multiples of sigma_effect
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C2 = scenario(
    "Multiple outliers: 3 of 20 from contaminating distribution",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 3L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C3 = scenario(
    "Heavy-tailed study effects: theta_i ~ t(df=3)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_t_effects",
      n_did      = 20L,
      effect_df  = 3
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C4 = scenario(
    "Heavy-tailed within-study errors: observations ~ t(df=5)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_t_likelihood",
      n_did      = 20L,
      within_df  = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C5 = scenario(
    "Outlier + small sample: single outlier among 5 studies",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 5L,
      n_outlier  = 1L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C6 = scenario(
    "Asymmetric contamination: outliers all biased in same direction",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 3L,
      outlier_shift = 5,
      outlier_direction = "positive"
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C7 = scenario(
    "Outlier in one design: 10 DiD + 10 PP, outliers only among PP",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_design_outliers",
      n_did      = 10L,
      n_pp       = 10L,
      n_outlier_pp = 2L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  C8 = scenario(
    "Heavy-tailed time trends: beta_i ~ t(df=3)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_t_trends",
      n_did      = 20L,
      trend_df   = 3
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  # ---------------------------------------------------------------------------
  # Category D: Assumption violations
  # ---------------------------------------------------------------------------

  D1 = scenario(
    "Design-specific effects ignored: RCT/PP shifted, design_effects=FALSE",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_design_offsets",
      n_did      = 10L,
      n_rct      = 10L,
      n_pp       = 10L,
      delta_rct  = 0.10,
      delta_pp   = -0.08
    )
  ),

  D2 = scenario(
    "Heterogeneous within-study variance: sigma_within drawn from Gamma",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_hetero_variance",
      n_did      = 20L,
      within_sd_shape = 4,     # Gamma shape
      within_sd_rate  = 4/0.12 # Gamma rate (mean = 0.12)
    )
  ),

  D3 = scenario(
    "Misspecified rho: true rho varies by study (Beta distribution)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_varying_rho",
      n_did      = 20L,
      rho_alpha  = 5,    # Beta(5, 5) => mean 0.5, moderate spread
      rho_beta   = 5
    )
  ),

  D4 = scenario(
    "Effect heterogeneity correlated with sample size",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_size_effect_corr",
      n_did      = 20L,
      size_effect_slope = 0.001 # larger studies have larger (less negative) effects
    )
  ),

  D5 = scenario(
    "Baseline imbalance in RCT with baseline_imbalance='fixed_zero'",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_rct_imbalance",
      n_did      = 10L,
      n_rct      = 10L,
      rct_gamma_mean = 0.05,   # mean baseline imbalance
      rct_gamma_sd   = 0.02
    ),
    fit = list(
      fn                 = "meta_did_general",
      baseline_imbalance = "fixed_zero"
    )
  ),

  # ---------------------------------------------------------------------------
  # Category E: Edge cases
  # ---------------------------------------------------------------------------

  E1 = scenario(
    "Extreme rho: pre-post correlation near 0 and near 1",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_extreme_rho",
      n_did      = 20L,
      rho_values = c(0.02, 0.98) # test both extremes
    )
  ),

  E2 = scenario(
    "Zero heterogeneity: all studies have identical true effect",
    dgp = list(sigma_effect = 0, sigma_trend = 0, baseline_sd = 0)
  ),

  E3 = scenario(
    "Large heterogeneity: tau >> typical effect size",
    dgp = list(sigma_effect = 0.30, sigma_trend = 0.10)
  ),

  E4 = scenario(
    "Unbalanced arms: n_treatment = 200, n_control = 30",
    dgp = list(n_control = 30L, n_treatment = 200L)
  ),

  E5 = scenario(
    "Single study per design: 1 DiD + 1 RCT + 1 PP",
    dgp = list(n_did = 1L, n_rct = 1L, n_pp = 1L)
  ),

  E6 = scenario(
    "Missing rho: hierarchical_rho imputes",
    dgp = list(n_did = 20L),
    fit = list(provide_rho = FALSE, hierarchical_rho = TRUE)
  ),

  # ---------------------------------------------------------------------------
  # Category A (continued): Sign and direction variants
  # All other parameters held at default; only effect/trend signs vary.
  # ---------------------------------------------------------------------------

  A14 = scenario(
    "Positive effect: true_effect = +0.15, trend = -0.04",
    dgp = list(true_effect = 0.15)
  ),

  A15 = scenario(
    "Null effect: true_effect = 0, trend = -0.04",
    dgp = list(true_effect = 0)
  ),

  A16 = scenario(
    "Positive trend: true_effect = -0.15, trend = +0.04",
    dgp = list(true_trend = 0.04)
  ),

  A17 = scenario(
    "Sign mismatch: negative effect (-0.15) with positive trend (+0.04), PP-heavy",
    dgp = list(n_did = 10L, n_pp = 20L, true_effect = -0.15, true_trend = 0.04)
  ),

  A18 = scenario(
    "Both positive: true_effect = +0.15, trend = +0.04",
    dgp = list(true_effect = 0.15, true_trend = 0.04)
  ),

  # ---------------------------------------------------------------------------
  # Category G: Bias source investigation
  #
  # G1-G5: Prior sensitivity (unnormalised, 20 DiD = same DGP as A9).
  #   Tests whether positive bias in unnormalised scenarios stems from
  #   prior shrinkage on treatment_effect_mean and time_trend_mean.
  #   Both priors widened together; if shrinkage is the cause, bias
  #   should diminish monotonically as prior sd increases.
  #
  # G6-G9: Jensen's inequality (normalised, 200 DiD — large N for sensitivity).
  #   Tests whether the negative bias in normalised scenarios arises from
  #   the ratio estimator: E[theta/baseline] != E[theta]/E[baseline].
  #   G6-G8 vary between-study baseline variation (baseline_sd).
  #   G9 tests within-study precision via a very large control group.
  # ---------------------------------------------------------------------------

  G1 = scenario(
    "Prior sensitivity: effect/trend prior sd = 1 (unnormalised, 20 DiD)",
    dgp = list(),
    fit = list(
      normalise = FALSE,
      priors = metadid::set_priors(
        treatment_effect_mean = metadid::normal(0, 1),
        time_trend_mean       = metadid::normal(0, 1)
      )
    )
  ),

  G2 = scenario(
    "Prior sensitivity: effect/trend prior sd = 5 (unnormalised, 20 DiD)",
    dgp = list(),
    fit = list(
      normalise = FALSE,
      priors = metadid::set_priors(
        treatment_effect_mean = metadid::normal(0, 5),
        time_trend_mean       = metadid::normal(0, 5)
      )
    )
  ),

  G3 = scenario(
    "Prior sensitivity: effect/trend prior sd = 10 (default, unnormalised, 20 DiD)",
    dgp = list(),
    fit = list(
      normalise = FALSE,
      priors = metadid::set_priors(
        treatment_effect_mean = metadid::normal(0, 10),
        time_trend_mean       = metadid::normal(0, 10)
      )
    )
  ),

  G4 = scenario(
    "Prior sensitivity: effect/trend prior sd = 50 (unnormalised, 20 DiD)",
    dgp = list(),
    fit = list(
      normalise = FALSE,
      priors = metadid::set_priors(
        treatment_effect_mean = metadid::normal(0, 50),
        time_trend_mean       = metadid::normal(0, 50)
      )
    )
  ),

  G5 = scenario(
    "Prior sensitivity: effect/trend prior sd = 100 (unnormalised, 20 DiD)",
    dgp = list(),
    fit = list(
      normalise = FALSE,
      priors = metadid::set_priors(
        treatment_effect_mean = metadid::normal(0, 100),
        time_trend_mean       = metadid::normal(0, 100)
      )
    )
  ),

  G6 = scenario(
    "Jensen's test: baseline_sd = 0, no between-study variation (200 DiD)",
    dgp = list(n_did = 200L, baseline_sd = 0)
  ),

  G7 = scenario(
    "Jensen's test: baseline_sd = 0.01, reduced variation (200 DiD)",
    dgp = list(n_did = 200L, baseline_sd = 0.01)
  ),

  G8 = scenario(
    "Jensen's test: baseline_sd = 0.10, amplified variation (100 DiD)",
    dgp = list(n_did = 100L, baseline_sd = 0.10)
  ),

  G9 = scenario(
    "Jensen's test: n_control = 5000, high within-study precision (100 DiD)",
    dgp = list(n_did = 100L, n_control = 5000L)
  ),

  # ---------------------------------------------------------------------------
  # Category H: Time trend distributional misspecification
  #
  # Each scenario uses 15 DiD + 15 PP studies (PP studies require trend
  # estimation for effect identification, making distributional shape matter).
  # All scenarios compare three models on the same simulated data:
  #   full_normal — meta_did with normal heterogeneity (default)
  #   full_robust — meta_did with Student-t heterogeneity
  #   naive       — meta_did_general with time_trend = "fixed_zero" and
  #                 baseline_imbalance = "fixed_zero"
  # ---------------------------------------------------------------------------

  H1 = scenario(
    "Skewed trends: beta_i ~ -LogNormal (mean = -0.04, log-SD = 1), 15 DiD + 15 PP",
    dgp = list(
      type             = "bespoke",
      bespoke_fn       = "simulate_lognormal_trends",
      n_did            = 15L,
      n_pp             = 15L,
      true_trend       = -0.04,   # E[beta_i] matches lognormal mean
      sigma_trend      = 0.052,   # approx SD of the log-normal
      lognormal_sigma  = 1.0      # log-scale SD; controls degree of skew
    ),
    compare = list(
      list(label = "full_normal", fn = "meta_did"),
      list(label = "full_robust", fn = "meta_did", robust_heterogeneity = TRUE),
      list(label = "naive",       fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  H2 = scenario(
    "Bimodal trends: 50/50 mixture N(-0.10, 0.01) + N(-0.01, 0.01), 15 DiD + 15 PP",
    dgp = list(
      type         = "bespoke",
      bespoke_fn   = "simulate_bimodal_trends",
      n_did        = 15L,
      n_pp         = 15L,
      true_trend   = -0.055,  # mixture mean
      sigma_trend  = 0.046,   # mixture SD
      trend_low    = -0.10,   # mean of low-trend component
      trend_high   = -0.01,   # mean of high-trend component
      mix_prob     = 0.5,     # probability of high-trend component
      sigma_within = 0.01     # within-component SD
    ),
    compare = list(
      list(label = "full_normal", fn = "meta_did"),
      list(label = "full_robust", fn = "meta_did", robust_heterogeneity = TRUE),
      list(label = "naive",       fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  H3 = scenario(
    "Large trend variance: sigma_trend = 0.10 (5x default), 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, sigma_trend = 0.10),
    compare = list(
      list(label = "full_normal", fn = "meta_did"),
      list(label = "full_robust", fn = "meta_did", robust_heterogeneity = TRUE),
      list(label = "naive",       fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  H4 = scenario(
    "Trend-size confounding: DiD n ~ U(30,70), PP n ~ U(150,250), 15 DiD + 15 PP",
    dgp = list(
      type             = "bespoke",
      bespoke_fn       = "simulate_trend_size_corr",
      n_did            = 15L,
      n_pp             = 15L,
      true_trend       = -0.04,   # trend at mean study size
      sigma_trend      = 0.02,    # residual trend SD after size adjustment
      trend_size_slope = 3e-4,    # beta units per unit of n_treatment above mean
      n_small_min      = 30L,     # DiD study size range
      n_small_max      = 70L,
      n_large_min      = 150L,    # PP study size range
      n_large_max      = 250L
    ),
    compare = list(
      list(label = "full_normal", fn = "meta_did"),
      list(label = "full_robust", fn = "meta_did", robust_heterogeneity = TRUE),
      list(label = "naive",       fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  )
)

# ---------------------------------------------------------------------------
# Helpers to extract scenario IDs by category
# ---------------------------------------------------------------------------

scenario_ids <- function(prefix) {
  ids <- grep(paste0("^", prefix, "\\d"), names(SCENARIO_CONFIGS), value = TRUE)
  sort(ids)
}

scenario_lookup <- function() {
  tibble::tibble(
    scenario_id = names(SCENARIO_CONFIGS),
    description = purrr::map_chr(SCENARIO_CONFIGS, "description")
  )
}

# Format a named list of overrides as a compact key=value string,
# omitting NULL values and collapsing to a single line.
.fmt_overrides <- function(overrides, defaults) {
  diffs <- overrides[!names(overrides) %in% c("type", "bespoke_fn", "covariates", "beta_cov")]
  diffs <- Filter(function(x) !is.null(x), diffs)
  # Keep only keys that differ from the defaults
  diffs <- diffs[purrr::map_lgl(names(diffs), function(k) {
    !identical(diffs[[k]], defaults[[k]])
  })]
  if (length(diffs) == 0) return("(defaults)")
  paste(names(diffs), purrr::map_chr(diffs, function(v) {
    if (is.numeric(v) && length(v) == 1) as.character(v)
    else if (is.character(v) && length(v) == 1) v
    else paste0("[", paste(v, collapse = ", "), "]")
  }), sep = " = ", collapse = "; ")
}

#' Scenario summary table for a given category
#'
#' Returns a data frame with one row per scenario showing the description,
#' DGP overrides from defaults, and fit overrides from defaults.
#'
#' @param category Character prefix, e.g. "A"
scenario_summary_table <- function(category) {
  ids <- stringr::str_sort(scenario_ids(category), numeric = TRUE)
  purrr::map_dfr(ids, function(id) {
    cfg <- SCENARIO_CONFIGS[[id]]
    dgp_str <- .fmt_overrides(cfg$dgp, default_dgp)
    fit_str <- .fmt_overrides(cfg$fit, default_fit)
    tibble::tibble(
      ID          = id,
      Description = cfg$description,
      `DGP overrides` = dgp_str,
      `Fit overrides` = fit_str
    )
  })
}

# ---------------------------------------------------------------------------
# Long-format scenario × model settings table
# ---------------------------------------------------------------------------

# Which programmatic comparators each category runs (see _targets.R).
# These mirror the matrix of tar_map_rep targets and must stay in sync.
.COMPARATOR_RULES <- list(
  A = c("naive", "robust"),
  B = character(),
  C = character(),
  D = c("naive", "robust"),
  E = c("naive", "robust"),
  F = c("robust"),
  G = character(),
  H = character()
)

# Coerce one config value (which may be NULL, a formula, a data.frame, a
# named list of priors, a vector, or a scalar) to a single string so it can
# sit in one CSV cell.
.flatten_value <- function(x) {
  if (is.null(x))                   return(NA_character_)
  if (inherits(x, "formula"))       return(paste(deparse(x), collapse = " "))
  if (is.data.frame(x))             return(paste0("data.frame(", paste(names(x), collapse = ","), ")"))
  if (is.list(x))                   return(paste(deparse(x), collapse = " "))
  if (length(x) == 0)               return(NA_character_)
  if (length(x) > 1)                return(paste(x, collapse = ";"))
  as.character(x)
}

# Flatten a named list (e.g. cfg$dgp or cfg$fit) into a named list of
# single-string values, prefixing each key. `defaults` is the corresponding
# defaults list so the column set is stable across scenarios.
.flatten_named_list <- function(values, prefix, defaults) {
  keys <- union(names(defaults), names(values))
  out  <- lapply(keys, function(k) .flatten_value(values[[k]]))
  names(out) <- paste0(prefix, "_", keys)
  out
}

# The fixed naive model config (mirrors fit_naive.R::run_naive_rep).
.build_naive_fit <- function(base_fit) {
  list(
    fn                    = "meta_did_general",
    normalise             = base_fit$normalise,
    baseline_latent_arm   = base_fit$baseline_latent_arm %||% "treatment",
    robust_heterogeneity  = FALSE,
    design_effects        = base_fit$design_effects,
    correlated_effects    = FALSE,
    hierarchical_rho      = base_fit$hierarchical_rho,
    time_trend            = "fixed_zero",
    baseline_imbalance    = "fixed_zero",
    pp_likelihood         = default_fit$pp_likelihood,
    covariates            = NULL,
    provide_rho           = default_fit$provide_rho,
    data_format           = "summary"
  )
}

.make_settings_row <- function(scenario_id, cfg, fit_for_row, model_label) {
  base <- list(
    scenario_id = scenario_id,
    category    = substr(scenario_id, 1L, 1L),
    description = cfg$description,
    model_label = model_label
  )
  dgp_flat <- .flatten_named_list(cfg$dgp, "dgp", default_dgp)
  fit_full <- modifyList(default_fit, fit_for_row)
  fit_flat <- .flatten_named_list(fit_full, "fit", default_fit)
  tibble::as_tibble(c(base, dgp_flat, fit_flat))
}

#' Build a long-format table with one row per (scenario_id, model) actually run
#'
#' Mirrors the model-execution logic in `_targets.R`. Every scenario contributes
#' one row per primary fit (the `compare` list, if present, expands into multiple
#' rows), plus a "naive" row for A/D/E scenarios whose base `data_format` is not
#' `"individual"`, plus a "robust" row for A/F/D/E scenarios whose base
#' `correlated_effects` is not `TRUE`. G scenarios contribute only their primary
#' fit (matching `run_g_rep`).
#'
#' Returned columns:
#'   scenario_id, category, description, model_label,
#'   dgp_*  — every key in `default_dgp` plus any scenario-specific keys,
#'   fit_*  — every key in `default_fit` plus any scenario-specific keys
#'           (e.g. `fit_priors` for G scenarios).
build_scenario_settings_table <- function() {
  rows <- list()
  for (sid in names(SCENARIO_CONFIGS)) {
    cfg         <- SCENARIO_CONFIGS[[sid]]
    category    <- substr(sid, 1L, 1L)
    comparators <- .COMPARATOR_RULES[[category]]
    if (is.null(comparators)) comparators <- character()

    # --- Primary fit(s) ---
    if (!is.null(cfg$compare)) {
      for (cmp in cfg$compare) {
        merged_fit <- modifyList(cfg$fit, cmp)
        label      <- if (!is.null(cmp$label)) cmp$label else "default"
        rows[[length(rows) + 1L]] <- .make_settings_row(sid, cfg, merged_fit, label)
      }
    } else {
      rows[[length(rows) + 1L]] <- .make_settings_row(sid, cfg, cfg$fit, "default")
    }

    # --- Naive comparator ---
    if ("naive" %in% comparators &&
        !isTRUE(cfg$fit$data_format == "individual")) {
      rows[[length(rows) + 1L]] <- .make_settings_row(
        sid, cfg, .build_naive_fit(cfg$fit), "naive"
      )
    }

    # --- Robust comparator ---
    if ("robust" %in% comparators &&
        !isTRUE(cfg$fit$correlated_effects)) {
      robust_fit <- modifyList(cfg$fit, list(robust_heterogeneity = TRUE))
      rows[[length(rows) + 1L]] <- .make_settings_row(sid, cfg, robust_fit, "robust")
    }
  }
  dplyr::bind_rows(rows)
}
