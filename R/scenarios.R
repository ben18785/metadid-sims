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
#     robust_heterogeneity, design_effects, correlated_effects
#     hierarchical_rho
#     time_trend, baseline_imbalance, pp_likelihood (for meta_did_general)
#     covariates      (formula or NULL)
#     provide_rho     Whether to include rho in summary data (default TRUE).
#                     TRUE = all studies report rho, FALSE = none, or a
#                     fraction p in [0,1] = the first round(p*n) studies report
#                     rho and the rest are left missing (for hierarchical
#                     imputation with an anchoring subset).
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
  fn                       = "meta_did",
  normalise                = TRUE,
  robust_heterogeneity     = FALSE,
  design_effects           = FALSE,
  correlated_effects       = FALSE,
  hierarchical_rho         = TRUE,
  time_trend               = "pooled",
  baseline_imbalance       = "estimated",
  pp_likelihood            = "differenced",
  covariates               = NULL,
  multiplicative_covariate = NULL,
  provide_rho              = TRUE,
  data_format              = "summary"
)

# Size of the DiD core that the composition sweeps build on. Exported so the
# panel code derives "studies added" instead of hard-coding the value.
K_CORE <- 5L    # category K (panel B): small core, default effect heterogeneity
S_CORE <- 15L   # category S (panel E): high effect heterogeneity

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
    # The design offsets must be injected into the data (dgp$delta_rct/pp,
    # read by simulate_design_offsets) AND declared as truths (true$delta_rct/pp,
    # read by build_true_params). The default `type = "metadid"` simulator has
    # no offset argument, so routing here is required for the RCT/PP arms to
    # actually carry the shift the fit is scored against.
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_design_offsets",
      n_did      = 35L, n_rct = 35L, n_pp = 30L,
      delta_rct  = 0.10, delta_pp = -0.08
    ),
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
    # Half the studies report rho; the hierarchical model imputes the rest by
    # borrowing strength from the reported ones. (Fitting with NO reported rho
    # is now an error in metadid — the hierarchy has nothing to anchor to — so
    # this scenario supplies an anchoring subset rather than none.)
    "Partial rho: hierarchical_rho imputes the missing half from the reported half",
    dgp = list(n_did = 20L),
    fit = list(provide_rho = 0.5, hierarchical_rho = TRUE)
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
  ),

  # ---------------------------------------------------------------------------
  # Category I: Multiplicative covariate scenarios
  #
  # Tests the multiplicative-covariate feature on the metadid
  # multiplicative-covars branch: an estimated per-level multiplier applied
  # to the population-mean linear predictor (treatment_effect_mean +
  # X_cov · beta_cov). Studies at the reference level keep the linear
  # predictor unchanged; studies at non-reference levels have it multiplied
  # by effect_multiplier[<level>], hierarchically distributed around the
  # multiplied mean.
  #
  # All I-category scenarios use the bespoke simulator
  # simulate_multiplicative_levels (see R/simulate.R) and pair modelled
  # (normalise = TRUE) vs raw (normalise = FALSE) fits via the compare
  # block. I5 and I6 additionally pair with-multiplier vs without-multiplier
  # fits (4 fits per rep) to quantify misspecification.
  # ---------------------------------------------------------------------------

  I1 = scenario(
    "Multiplicative covariate: binary balanced, 20 DiD",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 20L,
      level_assignments = rep_len(1:2, 20),
      level_multipliers = c(1, 0.6)
    ),
    fit = list(multiplicative_covariate = "level"),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  I2 = scenario(
    "Multiplicative covariate: three-level categorical, 30 DiD",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 30L,
      level_assignments = rep_len(1:3, 30),
      level_multipliers = c(1, 0.7, 0.4)
    ),
    fit = list(multiplicative_covariate = "level"),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  I3 = scenario(
    "Multiplicative + additive: binary multiplier × continuous dose, 30 DiD",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 30L,
      covariates        = data.frame(dose = seq(1, 4, length.out = 30)),
      beta_cov          = -0.04,
      level_assignments = rep_len(1:2, 30),
      level_multipliers = c(1, 0.6)
    ),
    fit = list(
      covariates               = ~ dose,
      multiplicative_covariate = "level"
    ),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  I4 = scenario(
    "Multiplicative covariate: mixed designs (10 DiD + 10 RCT + 10 PP), binary",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 10L,
      n_rct             = 10L,
      n_pp              = 10L,
      level_assignments = rep_len(1:2, 30),
      level_multipliers = c(1, 0.6)
    ),
    fit = list(multiplicative_covariate = "level"),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  # I5: truth has multiplicative structure; fit with vs without the multiplier
  # (each × modelled/raw). Quantifies the cost of omitting the feature.
  I5 = scenario(
    "Omit-multiplier misspecification: truth has binary multiplier, fit with vs without (× modelled/raw)",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 30L,
      level_assignments = rep_len(1:2, 30),
      level_multipliers = c(1, 0.6)
    ),
    compare = list(
      list(label = "with_modelled",    multiplicative_covariate = "level", normalise = TRUE),
      list(label = "with_raw",         multiplicative_covariate = "level", normalise = FALSE),
      list(label = "without_modelled", multiplicative_covariate = NULL,    normalise = TRUE),
      list(label = "without_raw",      multiplicative_covariate = NULL,    normalise = FALSE)
    )
  ),

  # I6: truth has NO multiplicative effect (all true multipliers = 1); fit
  # with a spurious multiplicative covariate vs without (each × modelled/raw).
  # Tests that adding the feature when unnecessary doesn't introduce bias.
  I6 = scenario(
    "Spurious multiplier: truth has NO multiplicative effect, fit with vs without (× modelled/raw)",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 30L,
      level_assignments = rep_len(1:2, 30),
      level_multipliers = c(1, 1)
    ),
    compare = list(
      list(label = "with_modelled",    multiplicative_covariate = "level", normalise = TRUE),
      list(label = "with_raw",         multiplicative_covariate = "level", normalise = FALSE),
      list(label = "without_modelled", multiplicative_covariate = NULL,    normalise = TRUE),
      list(label = "without_raw",      multiplicative_covariate = NULL,    normalise = FALSE)
    )
  ),

  I7 = scenario(
    "Multiplicative covariate, individual-level data, 20 DiD",
    dgp = list(
      type              = "bespoke",
      bespoke_fn        = "simulate_multiplicative_levels",
      n_did             = 20L,
      level_assignments = rep_len(1:2, 20),
      level_multipliers = c(1, 0.6)
    ),
    fit = list(
      data_format              = "individual",
      multiplicative_covariate = "level"
    ),
    # Only the modelled (normalised) arm. On the raw scale the individual-data +
    # multiplicative posterior is badly conditioned: the fit fails to converge
    # (Rhat ~2.1, ESS ~3) AND spends ~156 iterations/fit at max treedepth, which
    # made I7 ~3x slower than any other scenario for zero usable signal. The raw
    # multiplicative feature is already exercised by I1/I2/I8/I9's raw arms.
    compare = list(
      list(label = "modelled", normalise = TRUE)
    )
  ),

  # I8/I9: TWO multiplicative covariates fitted as a product via the formula
  # interface (~ level + level2). Each study's overall factor is the product
  # of its two per-covariate multipliers; assignments are crossed so every
  # cell of the (level × level2) table is occupied and the two factors are
  # separately identified.

  I8 = scenario(
    "Two-multiplier product: binary (1, 0.6) × binary (1, 1.4), 30 DiD",
    dgp = list(
      type               = "bespoke",
      bespoke_fn         = "simulate_multiplicative_levels",
      n_did              = 30L,
      level_assignments  = rep_len(1:2, 30),
      level_multipliers  = c(1, 0.6),
      level2_assignments = rep_len(rep(1:2, each = 2), 30),
      level2_multipliers = c(1, 1.4)
    ),
    fit = list(multiplicative_covariate = ~ level + level2),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  I9 = scenario(
    "Two-multiplier product: binary (1, 0.6) × three-level (1, 0.7, 0.4), 30 DiD",
    dgp = list(
      type               = "bespoke",
      bespoke_fn         = "simulate_multiplicative_levels",
      n_did              = 30L,
      level_assignments  = rep_len(1:2, 30),
      level_multipliers  = c(1, 0.6),
      level2_assignments = rep_len(rep(1:3, each = 2), 30),
      level2_multipliers = c(1, 0.7, 0.4)
    ),
    fit = list(multiplicative_covariate = ~ level + level2),
    compare = list(
      list(label = "modelled", normalise = TRUE),
      list(label = "raw",      normalise = FALSE)
    )
  ),

  # ---------------------------------------------------------------------------
  # Categories J-N: parameter sweeps for the paper's simulation figure.
  #
  #   J: trend-MEAN sweep (panel B)     - naive vs full bias as |mu_beta| grows
  #   K: composition sweep (panel D)    - posterior SD as studies are added
  #   L: contamination sweep (panel E)  - normal vs Student-t coverage
  #   M: trend-VARIANCE sweep (panel C) - mu_beta = 0; naive miscalibrates
  #   N: exchangeability sweep (panel F)- bias vs PP-DiD trend gap Delta
  #
  # These run at a reduced replication count (N_REPS_FIG in _targets.R)
  # while the figure design settles; raise N_REPS_FIG for the final run.
  # ---------------------------------------------------------------------------

  J1 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = 0, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = 0, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J2 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.025, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.025, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J3 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.05, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.05, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J4 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.075, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.075, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J5 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.1, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.1, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J6 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.125, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.125, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J7 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.15, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.15, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J8 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.175, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.175, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  J9 = scenario(
    "Figure panel B: trend-mean sweep, true_trend = -0.2, 5 DiD + 25 PP",
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = -0.2, sigma_trend = 0.01),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  # Core: 5 DiD studies. Panel B asks what adding lower-quality designs
  # buys when the DiD evidence base is SMALL, so the core is deliberately
  # thin; K_CORE is exported so the panel code derives "studies added"
  # rather than hard-coding it.
  K1 = scenario(
    "Figure panel B: 5 DiD + 0 RCT + 0 PP (0 incomplete studies added)",
    dgp = list(n_did = 5L, n_rct = 0L, n_pp = 0L)
  ),

  K2 = scenario(
    "Figure panel B: 5 DiD + 5 RCT + 5 PP (10 incomplete studies added)",
    dgp = list(n_did = 5L, n_rct = 5L, n_pp = 5L)
  ),

  K3 = scenario(
    "Figure panel B: 5 DiD + 10 RCT + 10 PP (20 incomplete studies added)",
    dgp = list(n_did = 5L, n_rct = 10L, n_pp = 10L)
  ),

  K4 = scenario(
    "Figure panel B: 5 DiD + 20 RCT + 20 PP (40 incomplete studies added)",
    dgp = list(n_did = 5L, n_rct = 20L, n_pp = 20L)
  ),

  K5 = scenario(
    "Figure panel B: 5 DiD + 40 RCT + 40 PP (80 incomplete studies added)",
    dgp = list(n_did = 5L, n_rct = 40L, n_pp = 40L)
  ),

  K6 = scenario(
    "Figure panel B reference: 15 DiD studies (10 DiD added to core of 5)",
    dgp = list(n_did = 15L)
  ),

  K7 = scenario(
    "Figure panel B reference: 25 DiD studies (20 DiD added to core of 5)",
    dgp = list(n_did = 25L)
  ),

  K8 = scenario(
    "Figure panel B reference: 45 DiD studies (40 DiD added to core of 5)",
    dgp = list(n_did = 45L)
  ),

  K9 = scenario(
    "Figure panel B reference: 85 DiD studies (80 DiD added to core of 5)",
    dgp = list(n_did = 85L)
  ),

  K10 = scenario(
    "Figure panel B split: 5 DiD + 10 RCT only (10 added)",
    dgp = list(n_did = 5L, n_rct = 10L, n_pp = 0L)
  ),

  K11 = scenario(
    "Figure panel B split: 5 DiD + 20 RCT only (20 added)",
    dgp = list(n_did = 5L, n_rct = 20L, n_pp = 0L)
  ),

  K12 = scenario(
    "Figure panel B split: 5 DiD + 40 RCT only (40 added)",
    dgp = list(n_did = 5L, n_rct = 40L, n_pp = 0L)
  ),

  K13 = scenario(
    "Figure panel B split: 5 DiD + 80 RCT only (80 added)",
    dgp = list(n_did = 5L, n_rct = 80L, n_pp = 0L)
  ),

  K14 = scenario(
    "Figure panel B split: 5 DiD + 10 PP only (10 added)",
    dgp = list(n_did = 5L, n_rct = 0L, n_pp = 10L)
  ),

  K15 = scenario(
    "Figure panel B split: 5 DiD + 20 PP only (20 added)",
    dgp = list(n_did = 5L, n_rct = 0L, n_pp = 20L)
  ),

  K16 = scenario(
    "Figure panel B split: 5 DiD + 40 PP only (40 added)",
    dgp = list(n_did = 5L, n_rct = 0L, n_pp = 40L)
  ),

  K17 = scenario(
    "Figure panel B split: 5 DiD + 80 PP only (80 added)",
    dgp = list(n_did = 5L, n_rct = 0L, n_pp = 80L)
  ),

  L1 = scenario(
    "Figure panel E: 0 of 20 studies outlying (5*sigma_effect shift)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 0L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  L2 = scenario(
    "Figure panel E: 1 of 20 studies outlying (5*sigma_effect shift)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 1L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  L3 = scenario(
    "Figure panel E: 3 of 20 studies outlying (5*sigma_effect shift)",
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

  L4 = scenario(
    "Figure panel E: 5 of 20 studies outlying (5*sigma_effect shift)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 5L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  L5 = scenario(
    "Figure panel E: 7 of 20 studies outlying (5*sigma_effect shift)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_with_outliers",
      n_did      = 20L,
      n_outlier  = 7L,
      outlier_shift = 5
    ),
    compare = list(
      list(label = "normal", fn = "meta_did", robust_heterogeneity = FALSE),
      list(label = "robust", fn = "meta_did", robust_heterogeneity = TRUE)
    )
  ),

  M1 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  M2 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0.02, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0.02),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  M3 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0.04, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0.04),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  M4 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0.06, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0.06),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  M5 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0.1, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0.1),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  M6 = scenario(
    "Figure panel C: zero-mean trends, sigma_trend = 0.14, 15 DiD + 15 PP",
    dgp = list(n_did = 15L, n_pp = 15L, true_trend = 0, sigma_trend = 0.14),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N1 = scenario(
    "Figure panel F: exchangeability gap Delta = 0 (DiD trend -0.04, PP trend -0.04)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.04,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N2 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.025 (DiD trend -0.04, PP trend -0.065)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.065,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N3 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.05 (DiD trend -0.04, PP trend -0.09)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.09,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N4 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.075 (DiD trend -0.04, PP trend -0.115)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.115,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N5 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.1 (DiD trend -0.04, PP trend -0.14)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.14,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N6 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.125 (DiD trend -0.04, PP trend -0.165)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.165,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  N7 = scenario(
    "Figure panel F: exchangeability gap Delta = 0.15 (DiD trend -0.04, PP trend -0.19)",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.19,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  # ---------------------------------------------------------------------------
  # Category O: trend crossover sweep (figure panel G).
  # DiD trend fixed at -0.04; PP trend swept THROUGH it. Full model unbiased
  # only at exchangeability (pp = did); naive unbiased only at pp = 0.
  # ---------------------------------------------------------------------------

  O1 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.2",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.2,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O2 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.15",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.15,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O3 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.1",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.1,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O4 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.07",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.07,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O5 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.04",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.04,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O6 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend -0.02",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = -0.02,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O7 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend 0",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = 0,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O8 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend 0.04",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = 0.04,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  O9 = scenario(
    "Figure panel G: crossover, DiD trend -0.04, PP trend 0.08",
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = -0.04, pp_trend = 0.08,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  ),

  # ---------------------------------------------------------------------------
  # Category P: RCT baseline-imbalance sweep (figure panel H).
  # The gamma analogue of the panel-B trend sweep: naive (gamma = 0) bias
  # grows at the RCT weight W_R; estimating imbalance with DiD borrowing
  # reduces (but does not remove) the slope - imbalance is only partially
  # identified from post-only data (cf. B5/B6).
  # ---------------------------------------------------------------------------

  P1 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P2 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.02",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.02,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P3 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.04",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.04,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P4 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.05",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.05,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P5 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.06",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.06,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P6 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.08",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.08,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
    )
  ),

  P7 = scenario(
    "Figure panel H: RCT baseline imbalance, gamma mean 0.1",
    dgp = list(
      type           = "bespoke",
      bespoke_fn     = "simulate_rct_imbalance",
      n_did          = 10L,
      n_rct          = 20L,
      rct_gamma_mean = 0.1,
      rct_gamma_sd   = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general", baseline_imbalance = "fixed_zero")
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

# ---------------------------------------------------------------------------
# Category Q: trend-plane grid (figure panels K and L)
#
# Every N and O scenario fixes did_trend = -0.04, so the whole published
# exchangeability evidence base is a single vertical line in the
# (did_trend, pp_trend) plane. Q fills the plane on a symmetric grid, which
# lets the geometry be measured rather than extrapolated: the full model
# should be unbiased along the DIAGONAL pp = did (it borrows the DiD trend, so
# it is right whenever that trend transfers, at any magnitude), and the naive
# model along the HORIZONTAL pp = 0. Any departure from that -- e.g. bias
# depending on the absolute trend level and not only on the difference -- is
# itself the finding.
# ---------------------------------------------------------------------------
Q_GRID <- expand.grid(
  did_trend = c(-0.12, -0.06, 0, 0.06, 0.12),
  pp_trend  = c(-0.12, -0.06, 0, 0.06, 0.12)
)

for (.i in seq_len(nrow(Q_GRID))) {
  .did <- Q_GRID$did_trend[.i]
  .pp  <- Q_GRID$pp_trend[.i]
  SCENARIO_CONFIGS[[paste0("Q", .i)]] <- scenario(
    sprintf("Figure panels K/L: trend plane, DiD trend %s, PP trend %s",
            .did, .pp),
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 10L, n_pp = 10L,
      did_trend = .did, pp_trend = .pp,
      sigma_trend = 0.02
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  )
}
rm(.i, .did, .pp)

# ---------------------------------------------------------------------------
# Category R: how many DiD studies buy back the premium (figure panel J)
#
# The pure null -- mu_beta = 0 AND sigma_beta = 0, so there is no trend at all
# and the naive model is exactly correctly specified. A grid over the number of
# DiD studies (which identify the trend) crossed with the number of PP studies
# (which need it borrowed) -- panel J plots the full/naive interval RATIO, so
# the growing total N cancels and what is left is the premium itself.
#
# NOTE this is the best case for the full model: with sigma_beta = 0 every DiD
# study reads the same trend, so it is learned as fast as possible. Under
# heterogeneous trends more DiD studies would be needed for the same premium.
# ---------------------------------------------------------------------------
R_GRID <- expand.grid(
  n_did = c(2L, 5L, 10L, 15L, 20L, 30L),
  n_pp  = c(5L, 20L, 50L)
)

for (.i in seq_len(nrow(R_GRID))) {
  .nd <- R_GRID$n_did[.i]
  .np <- R_GRID$n_pp[.i]
  SCENARIO_CONFIGS[[paste0("R", .i)]] <- scenario(
    sprintf("Figure panel J: pure null (no trend), %d DiD + %d PP", .nd, .np),
    dgp = list(n_did = .nd, n_pp = .np, true_trend = 0, sigma_trend = 0),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  )
}
rm(.i, .nd, .np)

# ---------------------------------------------------------------------------
# Category W: does the naive weight share depend on trend heterogeneity?
#
# The analytic bias surface (.q_weights in R/plots.R) gives the naive model an
# inverse-variance weight share built from 1 / (tau_theta^2 + s^2). tau_beta
# never enters it, so the supplementary heterogeneity map has a naive band of
# EXACTLY constant height down every column. That is an assumption, not a
# result: the naive model sees PP "effects" of theta_i + beta_i, whose between-
# study variance is tau_theta^2 + tau_beta^2, and it estimates its heterogeneity
# from the data rather than being told the true tau_theta.
#
# The expectation is that the error is small, because the naive model fits ONE
# common tau across both blocks (it cannot tell PP from DiD), so the inflation
# lands on numerator and denominator together and largely cancels in the share.
# But that is an expectation, and neither existing category can test it: M
# sweeps sigma_trend at true_trend = 0, where naive bias is ~0 and the share is
# unidentified, while Q varies the mean trends at sigma_trend fixed to 0.02.
#
# So: hold a divergent MEAN trend fixed and sweep tau_beta. With did_trend = 0
# the naive bias is W_naive * pp_trend / baseline, so the measured bias traces
# W_naive(tau_beta) up to a known constant. Flat confirms the closed form;
# declining means the fitted naive model downweights its contaminated PP block
# and the analytic surface overstates naive bias at large tau_beta.
#
# metadid rides along as the sharper test of the other half of the formula --
# its bias SHOULD fall with tau_beta, via the identification discount.
#
# RESULT (40 reps): naive bias falls 5.7% (se 1.9%, p = 0.02) across the sweep,
# so the closed form's exact invariance is wrong but only mildly. The mechanism
# is confirmed -- inflating BOTH blocks by tau_beta^2 (one common fitted tau)
# predicts -4.9%; inflating PP alone predicts -77%. Separately, measured bias
# is below prediction in both arms (naive 0.89x, metadid 0.73x), which is not
# yet explained; see `.q_weights` in R/plots.R.
# ---------------------------------------------------------------------------
W_TAU_BETA <- c(0.01, 0.02, 0.04, 0.06, 0.09)
W_PP_TREND <- 0.12   # matches the outer edge of the Q grid

for (.i in seq_along(W_TAU_BETA)) {
  .tb <- W_TAU_BETA[.i]
  SCENARIO_CONFIGS[[paste0("W", .i)]] <- scenario(
    sprintf("SI heterogeneity check: PP trend %s, sigma_trend %s",
            W_PP_TREND, .tb),
    dgp = list(
      type       = "bespoke",
      bespoke_fn = "simulate_divergent_trends",
      n_did = 5L, n_pp = 5L,
      did_trend = 0, pp_trend = W_PP_TREND,
      sigma_trend = .tb
    ),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  )
}
rm(.i, .tb)

# ---------------------------------------------------------------------------
# Category U: what each incomplete design is worth (figure panel E)
#
# Per-study sampling variance of the theta contribution, with n per arm,
# within-subject SD sigma and pre-post correlation rho:
#
#   DiD  (pre+post, both arms)   4 sigma^2 (1 - rho) / n
#   RCT  (post only, both arms)  2 sigma^2 / n
#   PP   (pre+post, treated)     2 sigma^2 (1 - rho) / n
#
# so, as information relative to one DiD study,
#
#   PP  / DiD = 2            (flat in rho)
#   RCT / DiD = 2 (1 - rho)  (crosses parity at rho = 0.5)
#
# The DEFAULT rho of 0.5 sits exactly on that crossover, which is why the
# existing K splits cannot rank DiD against post-only RCTs. Sweeping rho puts
# each design's exchange rate on a measurable axis; the gap between the
# measured curve and the analytic line above is the price of incompleteness
# (PP must borrow the trend, the RCT must accommodate baseline imbalance).
#
# Deliberately run at the DEFAULT sigma_effect: the exchange rate only means
# anything while sampling variance dominates. Under high effect heterogeneity
# every design converges to parity, which is what category S shows.
# ---------------------------------------------------------------------------
U_CORES <- c(5L, 20L)   # DiD anchor sizes: weak and strong
U_ADD   <- 20L          # studies added, of whichever design
# Three rho values, not five: the theory panels are closed-form curves that
# need no simulation, so the grid only has to support the theory-vs-empirical
# validation figure. 0.5 is kept because s2_DiD == s2_RCT exactly there, making
# the RCT parity point invariant to tau_theta; 0.2 and 0.8 bracket it.
U_RHO   <- c(0.2, 0.5, 0.8)

# Two anchor levels because PP's ceiling of 2 assumes the trend is KNOWN. In
# practice it is estimated from the anchor, and that uncertainty is charged to
# every PP study. With a weak anchor the measured PP curve should sit well
# below its analytic line; with a strong one it should approach it. That is
# what stops the panel reading as "just collect pre-post studies".
#
# The anchor size is encoded in the scenario id (U05_..., U20_...) rather than
# carried in the dgp, so nothing unexpected is passed through to the simulator.
for (.core in U_CORES) {
  for (.r in U_RHO) {
    .tag <- paste0(sprintf("%02d", .core), "_",
                   sub("\\.", "", format(.r, nsmall = 2)))
    .arms <- list(
      core = list(n_did = .core),
      did  = list(n_did = .core + U_ADD),
      rct  = list(n_did = .core, n_rct = U_ADD),
      pp   = list(n_did = .core, n_pp  = U_ADD)
    )
    for (.a in names(.arms)) {
      SCENARIO_CONFIGS[[paste0("U", .tag, "_", .a)]] <- scenario(
        sprintf("Figure panel E: anchor %d DiD, rho = %s, %s", .core, .r, .a),
        dgp = modifyList(list(rho = .r), .arms[[.a]])
      )
    }
  }
}
rm(.core, .r, .tag, .arms, .a)

# ---------------------------------------------------------------------------
# Category T: the M sweep at a PP-heavy composition (figure panel F)
#
# Identical to M (mu_beta = 0, sigma_beta swept) but with 5 DiD + 25 PP instead
# of 15 + 15. The full/naive RMSE ratio should separate EARLIER here: the more
# of the evidence base leans on PP studies, the more ignoring between-study
# trend variation costs. Same lever as panel C, so C and F tell one story.
# ---------------------------------------------------------------------------
T_SIGMA <- c(0, 0.02, 0.04, 0.06, 0.1, 0.14)

for (.i in seq_along(T_SIGMA)) {
  .sd <- T_SIGMA[.i]
  SCENARIO_CONFIGS[[paste0("T", .i)]] <- scenario(
    sprintf("Figure panel F: zero-mean trends, sigma_trend = %s, 5 DiD + 25 PP",
            .sd),
    dgp = list(n_did = 5L, n_pp = 25L, true_trend = 0, sigma_trend = .sd),
    compare = list(
      list(label = "full",  fn = "meta_did"),
      list(label = "naive", fn = "meta_did_general",
           time_trend = "fixed_zero", baseline_imbalance = "fixed_zero")
    )
  )
}
rm(.i, .sd)

# ---------------------------------------------------------------------------
# Category S: information from incomplete designs under HIGH effect
# heterogeneity (figure panel; K composition sweep at sigma_effect = 0.10)
#
# Theory (SI S3): tau_theta^2 inflates V_D and V_P equally, pushing the
# per-study exchange rate V_D/V_P toward 1 -- incomplete designs approach
# parity with DiD studies precisely when every study is a noisy reading of
# the effect. This repeats the K mixed/reference sweep at sigma_effect =
# 0.10 (default 0.03): the mixed curve should sit CLOSER to the DiD
# reference than in K, though everything is less precise absolutely.
# ---------------------------------------------------------------------------
S_ADDED <- c(0L, 10L, 20L, 40L, 80L)

for (.i in seq_along(S_ADDED)) {
  .k <- S_ADDED[.i]
  SCENARIO_CONFIGS[[paste0("S", .i)]] <- scenario(
    sprintf("Figure panel (high heterogeneity): %d DiD + %d RCT + %d PP, sigma_effect = 0.10",
            S_CORE, .k %/% 2L, .k %/% 2L),
    dgp = list(n_did = S_CORE, n_rct = .k %/% 2L, n_pp = .k %/% 2L,
               sigma_effect = 0.10)
  )
}
for (.i in seq_along(S_ADDED[-1])) {
  .k <- S_ADDED[-1][.i]
  SCENARIO_CONFIGS[[paste0("S", length(S_ADDED) + .i)]] <- scenario(
    sprintf("Figure panel (high heterogeneity) reference: %d DiD, sigma_effect = 0.10",
            S_CORE + .k),
    dgp = list(n_did = S_CORE + .k, sigma_effect = 0.10)
  )
}
rm(.i, .k)

# ---------------------------------------------------------------------------
# Category V: the K composition sweep at the S core size (figure panel B)
#
# K is (low tau_theta, small core) and S is (high tau_theta, large core), so on
# their own they confound effect heterogeneity with core size. V is the missing
# cell -- low tau_theta at the LARGE core -- which makes the comparison
# attributable:
#
#   K vs V   core size, at fixed low tau_theta
#   V vs S   effect heterogeneity, at fixed core of 15
#
# Same added-study grid as S so the three curves share an x axis.
# ---------------------------------------------------------------------------
V_CORE <- 15L

for (.i in seq_along(S_ADDED)) {
  .k <- S_ADDED[.i]
  SCENARIO_CONFIGS[[paste0("V", .i)]] <- scenario(
    sprintf("Figure panel B: %d DiD + %d RCT + %d PP, default sigma_effect",
            V_CORE, .k %/% 2L, .k %/% 2L),
    dgp = list(n_did = V_CORE, n_rct = .k %/% 2L, n_pp = .k %/% 2L)
  )
}
for (.i in seq_along(S_ADDED[-1])) {
  .k <- S_ADDED[-1][.i]
  SCENARIO_CONFIGS[[paste0("V", length(S_ADDED) + .i)]] <- scenario(
    sprintf("Figure panel B reference: %d DiD, default sigma_effect",
            V_CORE + .k),
    dgp = list(n_did = V_CORE + .k)
  )
}
rm(.i, .k)


# TRUE if a scenario fits any individual-level data — either directly
# (fit$data_format) or via a compare arm. Individual-data fits are markedly
# slower, so the pipeline runs them at a reduced replication count.
scenario_is_individual <- function(s) {
  cfg <- SCENARIO_CONFIGS[[s]]
  isTRUE(cfg$fit$data_format == "individual") ||
    (!is.null(cfg$compare) &&
       any(vapply(cfg$compare,
                  function(a) isTRUE(a$data_format == "individual"),
                  logical(1))))
}

scenario_lookup <- function() {
  tibble::tibble(
    scenario_id = names(SCENARIO_CONFIGS),
    description = purrr::map_chr(SCENARIO_CONFIGS, "description")
  )
}

# ---------------------------------------------------------------------------
# Expected deviations registry
#
# Some scenarios (or specific comparison arms / parameters) are DESIGNED to
# show bias, under-coverage, or convergence trouble: deliberate misspecification
# demonstrations, boundary/degenerate parameter values, or normalisation stress
# tests. Their deviations are not defects, so the executive summary separates
# them from genuine issues instead of burying real problems in the same list.
#
# A registry row matches a flagged result when `scenario_id` is equal and
# `model_label` / `parameter` are either NA (wildcard = any) or equal to the
# flagged row's value.
# ---------------------------------------------------------------------------

scenario_expectations <- function() {
  tibble::tribble(
    ~scenario_id, ~model_label,       ~parameter,             ~reason,
    # --- Deliberate misspecification: a wrong/naive model fitted on purpose ---
    "B1", "naive",            NA,                     "Naive arm zeros a large PP time trend — misspecified by design ('full should win')",
    "B3", "full",             NA,                     "Demonstrates pooled-trend failure when trends differ by design ('naive wins')",
    "B3", "naive",            NA,                     "Naive comparison arm — misspecified by design",
    # RCT baseline imbalance is not per-study identifiable (one post obs per arm);
    # the full model leans on the hierarchical prior informed by the few DiD
    # studies, so under-coverage here is the demonstrated identification limit.
    "B5", "full",             NA,                     "RCT baseline imbalance not per-study identifiable; full model relies on hierarchical borrowing from few DiD — under-coverage expected",
    "B5", "naive",            NA,                     "Naive arm ignores RCT baseline imbalance — biased by design",
    "B6", "full",             NA,                     "Shared baseline-difference hierarchy is dominated by the larger DiD imbalance and mis-applies it to RCT — bias expected by design",
    "B6", "naive",            NA,                     "Naive arm ignores baseline imbalance — biased by design",
    "H1", "naive",            NA,                     "Naive arm — misspecified by design",
    "H2", "naive",            NA,                     "Naive arm — misspecified by design",
    "H3", "naive",            NA,                     "Naive arm — misspecified by design",
    "H4", NA,                 NA,                     "Trend confounded with study size and design (DiD small, PP large); the pooled-trend model cannot separate them — all arms biased by design",
    "C1", "normal",           NA,                     "Normal heterogeneity fitted to outlier-contaminated truth — inflation expected (robust arm is the calibrated one)",
    "C2", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    "C3", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    "C4", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    "C5", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    # Asymmetric (one-directional) contamination shifts the population mean even
    # under robust heterogeneity, so BOTH arms are expected to fail here.
    "C6", NA,                 NA,                     "One-directional contamination shifts the mean even under robust heterogeneity — both arms biased by design",
    "C7", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    "C8", "normal",           NA,                     "Normal fit to heavy-tailed/contaminated truth — expected (see robust arm)",
    "D4", NA,                 NA,                     "Effect correlated with sample size — informative-sampling misspecification demo",
    "D5", NA,                 NA,                     "baseline_imbalance = fixed_zero on real imbalance — misspecification demo",
    "I5", "without_modelled", NA,                     "Multiplier omitted when truth has it — misspecification demo",
    "I5", "without_raw",      NA,                     "Multiplier omitted when truth has it — misspecification demo",
    # --- Boundary / degenerate parameter values ---
    "E2", NA,                 "treatment_effect_sd",  "Zero true heterogeneity (SD = 0 boundary) — a positive-constrained SD cannot cover it",
    "E2", NA,                 "time_trend_sd",        "Zero true heterogeneity (SD = 0 boundary) — a positive-constrained SD cannot cover it",
    "E5", NA,                 NA,                     "One study per design — between-study heterogeneity is unidentified by construction",
    # --- Deliberate normalisation (Jensen) stress tests ---
    "G8", NA,                 NA,                     "Jensen's test: large baseline_sd deliberately stresses the normalisation nonlinearity",
    "G9", NA,                 NA,                     "Jensen's test: high within-study precision deliberately stresses the normalisation nonlinearity"
  ) |>
    # --- Figure-sweep categories (J-N): arms misspecified by design ---
    dplyr::bind_rows(
      tibble::tibble(
        scenario_id = scenario_ids("J"), model_label = "naive",
        parameter = NA_character_,
        reason = "Naive arm zeros the swept time trend - misspecified by design (figure panel B)"
      ),
      tibble::tibble(
        scenario_id = scenario_ids("M"), model_label = "naive",
        parameter = NA_character_,
        reason = "Naive arm ignores swept zero-mean trend variability - miscalibration expected by design (figure panel C)"
      ),
      tibble::tibble(
        scenario_id = scenario_ids("L"), model_label = "normal",
        parameter = NA_character_,
        reason = "Normal heterogeneity fitted to outlier-contaminated truth - expected (figure panel E; robust arm is the calibrated one)"
      ),
      tibble::tibble(
        scenario_id = scenario_ids("N"), model_label = NA_character_,
        parameter = NA_character_,
        reason = "Trend exchangeability deliberately violated - bias expected in both arms, bounded for the full model (figure panel F)"
      ),
      tibble::tibble(
        scenario_id = scenario_ids("O"), model_label = NA_character_,
        parameter = NA_character_,
        reason = "Crossover sweep deliberately violates trend exchangeability on both sides of the matched point (figure panel G)"
      ),
      tibble::tibble(
        scenario_id = scenario_ids("P"), model_label = NA_character_,
        parameter = NA_character_,
        reason = "RCT imbalance sweep: naive arm misspecified by design; full arm only partially identifies imbalance (figure panel H)"
      ),
      tibble::tibble(
        scenario_id = "M1", model_label = NA_character_,
        parameter = c("time_trend_sd"),
        reason = "Zero true trend heterogeneity (SD = 0 boundary) - a positive-constrained SD cannot cover it"
      )
    )
}

# Annotate a data frame of flagged results (must contain scenario_id,
# model_label, parameter) with `expected` (logical) and `reason` (chr),
# resolved against scenario_expectations() with NA-as-wildcard matching.
tag_expectations <- function(flagged) {
  reg <- scenario_expectations()
  hits <- purrr::pmap(
    list(flagged$scenario_id, flagged$model_label, flagged$parameter),
    function(sid, lab, par) {
      m <- reg[reg$scenario_id == sid &
                 (is.na(reg$model_label) | reg$model_label == lab) &
                 (is.na(reg$parameter)   | reg$parameter   == par), ]
      if (nrow(m) == 0) c(expected = FALSE, reason = NA_character_)
      else c(expected = TRUE, reason = m$reason[[1]])
    }
  )
  flagged$expected <- vapply(hits, function(h) as.logical(h[["expected"]]), logical(1))
  flagged$reason   <- vapply(hits, function(h) h[["reason"]], character(1))
  flagged
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
  H = character(),
  # Category I scenarios manage their own modelled/raw and with/without
  # comparators via per-scenario `compare` blocks, so no programmatic
  # cross-cutting comparators are added.
  I = character()
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
