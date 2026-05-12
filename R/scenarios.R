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
#     normalise_by_baseline
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
  normalise_by_baseline = TRUE,
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
    "Very large sample: DiD-only, 200 studies",
    dgp = list(n_did = 200L)
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
    "Correlated effects: 30 DiD, rho_effect_trend = 0.6",
    dgp = list(n_did = 30L, rho_effect_trend = 0.6, baseline_sd = 0),
    fit = list(fn = "meta_did_general", correlated_effects = TRUE)
  ),

  A9 = scenario(
    "Unnormalised DiD-only, 20 studies, raw scale",
    fit = list(normalise_by_baseline = FALSE)
  ),

  A10 = scenario(
    "Unnormalised mixed: 10 DiD + 10 RCT + 10 PP, raw scale",
    dgp = list(n_did = 10L, n_rct = 10L, n_pp = 10L),
    fit = list(normalise_by_baseline = FALSE)
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
    "Large mixed, normalised, full model: 70 DiD + 70 RCT + 60 PP",
    dgp = list(n_did = 70L, n_rct = 70L, n_pp = 60L)
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
    fit = list(normalise_by_baseline = FALSE)
  ),

  F4 = scenario(
    "Large DiD + RCT only: 100 DiD + 100 RCT",
    dgp = list(n_did = 100L, n_rct = 100L)
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
    "Large mixed, design effects: 70 DiD + 70 RCT + 60 PP with offsets",
    dgp = list(n_did = 70L, n_rct = 70L, n_pp = 60L),
    fit = list(design_effects = TRUE),
    true = list(delta_rct = 0.10, delta_pp = -0.08)
  ),

  F8 = scenario(
    "Large DiD-only, correlated effects: 200 studies, rho_effect_trend = 0.6",
    dgp = list(n_did = 200L, rho_effect_trend = 0.6, baseline_sd = 0),
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
    "Correlated vs independent effects: rho_effect_trend = 0.6",
    dgp = list(n_did = 30L, rho_effect_trend = 0.6, baseline_sd = 0),
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
