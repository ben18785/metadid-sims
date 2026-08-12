# Summary visualisation functions for the metadid validation study.

library(ggplot2)
library(dplyr)
library(gridExtra)
library(kableExtra)

# ===========================================================================
# Coverage calibration plot
# ===========================================================================

#' Coverage plot: empirical vs nominal 90% coverage by scenario and parameter
#'
#' @param agg_results Aggregated results from aggregate_scenario()
#' @param category Optional character prefix to filter scenarios (e.g., "A")
plot_coverage <- function(agg_results, category = NULL) {
  df <- agg_results
  if (!is.null(category)) {
    df <- df |> filter(grepl(paste0("^", category), scenario_id))
  }

  ggplot(df, aes(x = scenario_id, y = empirical_coverage, fill = parameter)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey40") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(
      x = "Scenario",
      y = "Empirical coverage (90% CI)",
      fill = "Parameter"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ===========================================================================
# Bias distribution plot
# ===========================================================================

#' Bias distribution across replications, faceted by parameter
#'
#' @param rep_results Per-replication results (stacked assess_one outputs)
#' @param category Optional character prefix to filter scenarios
plot_bias <- function(rep_results, category = NULL) {
  df <- rep_results
  if (!is.null(category)) {
    df <- df |> filter(grepl(paste0("^", category), scenario_id))
  }

  ggplot(df, aes(x = scenario_id, y = bias)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(x = "Scenario", y = "Bias (posterior mean - truth)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ===========================================================================
# Model comparison plot (for B-category)
# ===========================================================================

#' Compare bias between model configurations on the same data
#'
#' @param rep_results Per-replication results with model_label column
#' @param param Parameter name to compare (default: "treatment_effect_mean")
plot_comparison <- function(rep_results, param = "treatment_effect_mean") {
  df <- rep_results |>
    filter(parameter == param)

  ggplot(df, aes(x = model_label, y = abs(bias), fill = model_label)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ scenario_id, scales = "free_y") +
    labs(
      x = "Model configuration",
      y = paste0("|Bias| for ", param),
      fill = "Model"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ===========================================================================
# RMSE comparison across scenarios
# ===========================================================================

#' RMSE by scenario for a given parameter
#'
#' @param agg_results Aggregated results from aggregate_scenario()
#' @param param Parameter name (default: "treatment_effect_mean")
plot_rmse <- function(agg_results, param = "treatment_effect_mean") {
  df <- agg_results |>
    filter(parameter == param)

  ggplot(df, aes(x = scenario_id, y = rmse, fill = model_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    labs(
      x = "Scenario",
      y = paste0("RMSE for ", param),
      fill = "Model"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ===========================================================================
# Convergence diagnostic plot
# ===========================================================================

# ===========================================================================
# Three-panel diagnostic plot (% bias, coverage, Rhat)
# ===========================================================================

#' Three-panel diagnostic plot: % bias, coverage, and Rhat across replicates
#'
#' @param rep_results  Per-replication results (stacked assess_one outputs)
#' @param agg_results  Aggregated results from aggregate_scenario()
#' @param lookup       Scenario lookup table from scenario_lookup()
#' @param category     Character prefix to filter scenarios (e.g., "F")
#' @param param        Parameter to plot (default: "treatment_effect_mean")
plot_diagnostics <- function(rep_results, agg_results, lookup,
                             category, param = "treatment_effect_mean",
                             exclude = NULL) {

  cat_pattern <- paste0("^", category, "\\d")

  label_data <- lookup |>
    filter(grepl(cat_pattern, scenario_id)) |>
    mutate(
      label = paste0(scenario_id, ": ", stringr::str_wrap(description, 35)),
      # Natural numeric ordering: A1, A2, ..., A10, A11
      label = forcats::fct_reorder(
        label,
        stringr::str_extract(scenario_id, "\\d+") |> as.integer(),
        .desc = TRUE
      )
    )

  # Normalise model labels per scenario: single-model scenarios get NA so they
  # plot without colour; comparison scenarios retain their named labels.
  normalise_labels <- function(df) {
    df |>
      group_by(scenario_id) |>
      mutate(model_label = if (n_distinct(model_label) > 1) model_label
                           else if (model_label[[1]] == "default") "normal"
                           else model_label) |>
      ungroup()
  }

  rep_df <- rep_results |>
    filter(parameter == param, grepl(cat_pattern, scenario_id),
           !scenario_id %in% exclude) |>
    left_join(label_data, by = "scenario_id") |>
    normalise_labels() |>
    mutate(pct_bias = dplyr::if_else(true_value != 0,
                                     100 * bias / abs(true_value),
                                     bias))

  agg_df <- agg_results |>
    filter(parameter == param, grepl(cat_pattern, scenario_id),
           !scenario_id %in% exclude) |>
    left_join(label_data, by = "scenario_id") |>
    normalise_labels() |>
    mutate(
      successes = round(empirical_coverage * n_reps),
      ci_lo     = qbeta(0.025, successes,     n_reps - successes + 1),
      ci_hi     = qbeta(0.975, successes + 1, n_reps - successes)
    )

  no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                axis.title.y = element_blank())

  p1 <- ggplot(rep_df, aes(x = pct_bias, y = label, colour = model_label)) +
    geom_boxplot(width = 0.5, position = position_dodge(0.6)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = "Bias (%; raw bias where true = 0)", y = NULL, colour = "Model")

  p2 <- ggplot(agg_df, aes(x = empirical_coverage, y = label, colour = model_label)) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(0.6)) +
    geom_point(position = position_dodge(0.6)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Coverage", colour = "Model") +
    no_y

  p3 <- ggplot(rep_df, aes(x = max_rhat, y = label, colour = model_label)) +
    geom_boxplot(width = 0.5, position = position_dodge(0.6)) +
    geom_vline(xintercept = 1.05, linetype = "dashed", colour = "red") +
    labs(x = "Max Rhat", colour = "Model") +
    no_y

  grid.arrange(p1, p2, p3, nrow = 1, widths = c(3, 1.5, 3),
               top = paste0(category, " scenarios: ", param, " diagnostics"))
}

# ===========================================================================
# Results table with collapsed scenario descriptions
# ===========================================================================

#' Aggregated results table with one description block per scenario
#'
#' @param agg_results  Aggregated results from aggregate_scenario()
#' @param lookup       Scenario lookup table from scenario_lookup()
#' @param category     Character prefix to filter scenarios (e.g., "A")
#' @param caption      Table caption string
results_table <- function(agg_results, lookup, category, caption = NULL) {
  if (is.null(caption))
    caption <- paste0("Category ", category, ": results")

  # Summarise key model assumptions as a readable string
  fmt_assumptions <- function(fit) {
    `%||%` <- rlang::`%||%`
    parts <- c(
      if (isTRUE(fit$normalise))  "normalised"       else "unnormalised",
      if (isTRUE(fit$robust_heterogeneity))   "robust heterogeneity" else "normal heterogeneity",
      if (isTRUE(fit$design_effects))         "design effects"   else NULL,
      if (isTRUE(fit$correlated_effects))     "correlated effects" else NULL,
      switch(fit$time_trend %||% "pooled",
             pooled     = "pooled trend",
             fixed_zero = "zero trend",
             fit$time_trend),
      switch(fit$baseline_imbalance %||% "estimated",
             estimated  = "estimated imbalance",
             fixed_zero = "zero imbalance",
             fit$baseline_imbalance),
      if (!isTRUE(fit$hierarchical_rho))      "fixed rho"        else NULL,
      if (fit$data_format == "individual")    "individual data"  else NULL
    )
    paste(parts, collapse = "; ")
  }

  model_info <- purrr::map_dfr(
    stringr::str_sort(
      grep(paste0("^", category, "\\d"), names(SCENARIO_CONFIGS), value = TRUE),
      numeric = TRUE
    ),
    function(id) {
      tibble::tibble(
        scenario_id  = id,
        model        = fmt_assumptions(SCENARIO_CONFIGS[[id]]$fit)
      )
    }
  )

  df <- agg_results |>
    filter(grepl(paste0("^", category, "\\d"), scenario_id)) |>
    left_join(lookup,     by = "scenario_id") |>
    left_join(model_info, by = "scenario_id") |>
    # Natural sort: A1, A2, ..., A10, A11
    mutate(scenario_id = forcats::fct_reorder(
      scenario_id, stringr::str_extract(scenario_id, "\\d+") |> as.integer()
    )) |>
    arrange(scenario_id, model_label, parameter) |>
    select(
      `Scenario`    = scenario_id,
      `Description` = description,
      `Model`       = model,
      `Fit`         = model_label,
      `Parameter`   = parameter,
      `Coverage`    = empirical_coverage,
      `Bias`        = mean_bias,
      `RMSE`        = rmse
    )

  df |>
    kbl(digits = 3, caption = caption) |>
    collapse_rows(columns = 1:3, valign = "top") |>
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
}

#' Proportion of fits with Rhat > 1.05, by scenario
#'
#' @param agg_results Aggregated results from aggregate_scenario()
plot_convergence <- function(agg_results) {
  df <- agg_results |>
    filter(parameter == "treatment_effect_mean") |>
    select(scenario_id, model_label, pct_rhat_bad) |>
    distinct()

  ggplot(df, aes(x = scenario_id, y = pct_rhat_bad, fill = model_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = "Scenario",
      y = "% fits with Rhat > 1.05",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Sampler-divergence diagnostics by scenario
#'
#' Two stacked panels: (top) mean number of divergent transitions per fit,
#' (bottom) percentage of reps with any divergent transitions. Both are
#' useful — the mean indicates magnitude when divergences do occur, the
#' percentage indicates how widespread the problem is across replicate fits.
#'
#' @param agg_results Aggregated results from aggregate_scenario()
plot_divergences <- function(agg_results) {
  df <- agg_results |>
    filter(parameter == "treatment_effect_mean") |>
    select(scenario_id, model_label,
           mean_n_divergent, pct_reps_with_divergence,
           pct_reps_many_divergence) |>
    distinct()

  p_mean <- ggplot(df, aes(x = scenario_id, y = mean_n_divergent, fill = model_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    labs(x = NULL, y = "Mean divergent transitions per fit", fill = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_pct <- ggplot(df, aes(x = scenario_id, y = pct_reps_with_divergence, fill = model_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Scenario", y = "% reps with > 0 divergences", fill = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  gridExtra::grid.arrange(p_mean, p_pct, nrow = 2,
                           top = "Sampler divergences by scenario")
}


# ===========================================================================
# Paper figure (categories A/F calibration + J-P sweeps)
# ===========================================================================
#
# Builds the multi-panel simulation figure for the paper from aggregated
# results. Sweep x-values are derived from SCENARIO_CONFIGS so the figure
# stays correct if the grids change. All effect-scale quantities are on the
# baseline-normalised scale (raw sweep values divided by the DGP
# baseline_mean).

.fig_base <- function(sid) {
  b <- SCENARIO_CONFIGS[[sid]]$dgp$baseline_mean
  if (is.null(b)) 0.45 else b
}

.fig_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 10)
    )
}

.fig_model_cols <- c(
  full        = "#0072B2", naive = "#D55E00",
  normal      = "#0072B2", robust = "#009E73",
  # `homog` is the middle rung of the M-category nested ladder (mu_beta free,
  # sigma_beta pinned at ~0). Shares a hue with `robust`, but the two never
  # appear in the same panel.
  homog       = "#009E73"
)

# Diverging fill for signed bias surfaces (panels K/L). Deliberately NOT the
# categorical model hues: here colour encodes sign and magnitude, not identity.
# Default tolerable |bias| for the panel-I model-choice map, as a fraction of
# the true effect. A judgement call, not a property of the data: 0.10 condemned
# most of the trend plane as unusable, and the map is flat from 0.05 to 0.15,
# so 0.15 is the least strict value that still discriminates. Defined once so
# the entry points below cannot drift apart; override per call where needed.
Q_TOL_FRAC <- 0.15

# TRUE when a pivoted frame actually has both model arms. Panels are built
# eagerly, so without this a category absent from `all_agg` (e.g. previewing a
# partial figure_panels.csv locally) would error out the whole panel list.
.has_arms <- function(df) {
  all(c("full", "naive") %in% names(df)) && nrow(df) > 0
}

.fig_bias_fill <- function(limit) {
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-limit, limit),
    name = expression("Bias in " * mu[theta])
  )
}

#' Multi-panel paper figure from aggregated results
#'
# ---------------------------------------------------------------------------
# Category U: per-design information helpers
#
# Shared by the theory panels (Figure 1 panel E, supplementary) and the
# theory-vs-empirical validation figures, so both read one implementation.
# Equations referenced below are from docs/design-information-derivation.docx.
# ---------------------------------------------------------------------------

# anchor size and arm live in the scenario id (U05_050_pp), not the dgp
.u_meta <- function(ids) {
  purrr::map_dfr(ids, function(id) {
    tibble::tibble(
      scenario_id = id,
      anchor = as.integer(sub("^U([0-9]+)_.*$", "\\1", id)),
      arm    = sub("^.*_", "", id),
      rho    = SCENARIO_CONFIGS[[id]]$dgp$rho)
  })
}

# equation (6.5): (tau-aware ratio) x (identification discount)
.u_theory <- function(ids, anchors, rho = seq(0.15, 0.85, length.out = 80)) {
  if (!length(ids)) return(tibble::tibble())
  g <- SCENARIO_CONFIGS[[ids[1]]]$dgp
  sig <- g$within_sd; nn <- g$n_treatment
  tau <- g$sigma_effect; tau_b <- g$sigma_trend; k <- U_ADD
  tidyr::expand_grid(rho = rho, design = c("RCT (post-only)", "PP"),
                     anchor = anchors) |>
    dplyr::mutate(
      s2D  = 4 * sig^2 * (1 - rho) / nn,
      s2X  = ifelse(design == "PP", 2 * sig^2 * (1 - rho) / nn, 2 * sig^2 / nn),
      ID   = 1 / (tau^2 + s2D), IX = 1 / (tau^2 + s2X),
      i_nu = ifelse(design == "PP", 1 / (tau_b^2 + 2 * sig^2 * (1 - rho) / nn),
                                    1 / (2 * sig^2 / nn)),
      ratio = (IX / ID) * (anchor * i_nu) / (anchor * i_nu + k * IX))
}

# equation (5.2): the ceiling, i.e. the nuisance treated as known
.u_ceiling <- function(ids, rho = seq(0.15, 0.85, length.out = 80)) {
  if (!length(ids)) return(tibble::tibble())
  g <- SCENARIO_CONFIGS[[ids[1]]]$dgp
  sig <- g$within_sd; nn <- g$n_treatment; tau <- g$sigma_effect
  tidyr::expand_grid(rho = rho, design = c("RCT (post-only)", "PP")) |>
    dplyr::mutate(
      s2D = 4 * sig^2 * (1 - rho) / nn,
      s2X = ifelse(design == "PP", 2 * sig^2 * (1 - rho) / nn, 2 * sig^2 / nn),
      ratio = (tau^2 + s2D) / (tau^2 + s2X))
}

# measured exchange rate: ratio of precision INCREMENTS, so the prior cancels
.u_measured <- function(te, ids) {
  empty <- tibble::tibble(rho = numeric(), anchor = numeric(),
                          design = character(), ratio = numeric())
  if (!length(ids)) return(empty)
  d <- te |>
    dplyr::filter(scenario_id %in% ids) |>
    dplyr::inner_join(.u_meta(ids), by = "scenario_id")
  # a non-converged `core` arm sits in the denominator of every ratio at its
  # grid point, so drop it and leave the point missing rather than wrong
  if ("max_rhat_worst" %in% names(d)) {
    d <- dplyr::filter(d, is.na(max_rhat_worst) | max_rhat_worst <= 1.05)
  }
  w <- d |>
    dplyr::mutate(prec = 1 / mean_posterior_sd^2) |>
    dplyr::select(anchor, rho, arm, prec) |>
    tidyr::pivot_wider(names_from = arm, values_from = prec)
  if (!all(c("core", "did", "rct", "pp") %in% names(w)) || !nrow(w)) return(empty)
  w |>
    dplyr::transmute(rho, anchor,
                     `RCT (post-only)` = (rct - core) / (did - core),
                     `PP`              = (pp  - core) / (did - core)) |>
    tidyr::pivot_longer(c(-rho, -anchor), names_to = "design",
                        values_to = "ratio")
}

#' @param all_agg Aggregated results tibble (scenario x model x parameter).
#' @return Named list of ggplot objects, one per panel letter.
#'
#' Split out from plot_paper_figure() so a single panel can be rebuilt and
#' previewed locally from a downloaded figure_panels.csv, without rerunning
#' the pipeline:  paper_panels(readr::read_csv("figure_panels.csv"))$I
#' Panels whose category is absent from `all_agg` come back empty rather than
#' erroring, so a partial CSV still lets you inspect the panels it does cover.
paper_panels <- function(all_agg, paired = NULL, q_tol_frac = Q_TOL_FRAC) {
  te <- all_agg |> dplyr::filter(parameter == "treatment_effect_mean")

  # --- A: calibration coverage across A and F scenarios -------------------
  pa <- te |>
    dplyr::filter(grepl("^[AF][0-9]", scenario_id)) |>
    dplyr::mutate(
      se = sqrt(pmax(empirical_coverage * (1 - empirical_coverage), 1e-9) / n_reps),
      lo = pmax(empirical_coverage - 1.96 * se, 0),
      hi = pmin(empirical_coverage + 1.96 * se, 1),
      scenario_id = factor(
        scenario_id,
        levels = stringr::str_sort(unique(scenario_id), numeric = TRUE)
      )
    )
  # A/F rows are absent from the figure workflow (calibration is a nightly
  # concern), so guard the empty case rather than warning on every build.
  band <- if (nrow(pa)) {
    1.96 * sqrt(0.9 * 0.1 / max(pa$n_reps, na.rm = TRUE))
  } else NA_real_
  gA <- ggplot2::ggplot(pa, ggplot2::aes(scenario_id, empirical_coverage)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = 0.9 - band, ymax = min(0.9 + band, 1),
                      alpha = 0.15) +
    ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi),
                             colour = "#0072B2", size = 0.2) +
    ggplot2::labs(title = "A  Calibration (A/F scenarios)",
                  x = NULL, y = "Coverage of 90% CrI") +
    .fig_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 5.5,
                                                       vjust = 0.5))

  # --- helper: sweep panel of a summary statistic vs a DGP-derived x ------
  sweep_data <- function(prefix, xfun) {
    ids <- scenario_ids(prefix)
    xv <- purrr::map_dbl(ids, function(id) xfun(SCENARIO_CONFIGS[[id]]$dgp, id))
    te |>
      dplyr::filter(scenario_id %in% ids) |>
      dplyr::left_join(tibble::tibble(scenario_id = ids, x = xv),
                       by = "scenario_id")
  }

  # --- B: trend-mean sweep (J) --------------------------------------------
  pb <- sweep_data("J", function(d, id) abs(d$true_trend) / .fig_base(id))
  gB <- ggplot2::ggplot(pb, ggplot2::aes(x, mean_bias, colour = model_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "B  Systematic trends: bias (J)",
                  x = expression("True trend magnitude " * "|" * mu[beta] * "|" *
                                   " (normalised)"),
                  y = expression("Bias in " * mu[theta])) +
    .fig_theme()

  # --- C: trend-variability sweep at zero mean (M) ------------------------
  pc <- sweep_data("M", function(d, id) d$sigma_trend / 0.02) |>
    dplyr::mutate(
      se = sqrt(pmax(empirical_coverage * (1 - empirical_coverage), 1e-9) / n_reps),
      lo = pmax(empirical_coverage - 1.96 * se, 0),
      hi = pmin(empirical_coverage + 1.96 * se, 1)
    )
  gC <- ggplot2::ggplot(pc, ggplot2::aes(x, empirical_coverage,
                                         colour = model_label)) +
    ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi), size = 0.25,
                             position = ggplot2::position_dodge(width = 0.25)) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "C  Variable trends, zero mean: coverage (M)",
                  x = expression(tau[beta] * " (multiples of default), " *
                                   mu[beta] * " = 0"),
                  y = "Coverage of 90% CrI") +
    .fig_theme()

  # --- D: composition sweeps at low and high effect heterogeneity ---------
  # K and S combined. Colour separates the two heterogeneity settings, the
  # linetype keeps the within-setting contrast (mixed vs the DiD reference).
  # The message is the GAP between a solid line and the dashed line of the
  # same colour: wide at low tau_theta (an incomplete study is worth less
  # than a DiD study), closed at high tau_theta (parity).
  .composition <- function(prefix, core, setting) {
    ids <- scenario_ids(prefix)
    info <- purrr::map_dfr(ids, function(id) {
      d <- SCENARIO_CONFIGS[[id]]$dgp
      tibble::tibble(
        scenario_id = id,
        added = (d$n_did - core) + d$n_rct + d$n_pp,
        curve = dplyr::case_when(
          d$n_rct > 0 & d$n_pp > 0 ~ "RCT + PP",
          d$n_rct > 0              ~ "RCT only",
          d$n_pp  > 0              ~ "PP only",
          d$n_did > core           ~ "DiD (reference)",
          TRUE                     ~ "core"
        )
      )
    })
    info <- dplyr::filter(info, curve %in% c("RCT + PP", "DiD (reference)",
                                             "core"))
    keep <- setdiff(unique(info$curve), "core")
    all <- dplyr::bind_rows(
      dplyr::filter(info, curve != "core"),
      tidyr::crossing(dplyr::select(dplyr::filter(info, curve == "core"),
                                    -curve), curve = keep)
    )
    te |>
      dplyr::filter(scenario_id %in% ids) |>
      dplyr::inner_join(all, by = "scenario_id",
                        relationship = "many-to-many") |>
      dplyr::mutate(setting = setting)
  }
  # Both settings differ in TWO respects -- effect heterogeneity and core size
  # -- so the legend names both explicitly rather than saying "low"/"high".
  k_se <- SCENARIO_CONFIGS[[scenario_ids("K")[1]]]$dgp$sigma_effect
  s_se <- SCENARIO_CONFIGS[[scenario_ids("S")[1]]]$dgp$sigma_effect
  pd <- dplyr::bind_rows(
    .composition("K", K_CORE, "low"),
    .composition("S", S_CORE, "high")
  ) |>
    dplyr::mutate(setting = factor(setting, levels = c("low", "high")))
  pd_labs <- list(
    bquote(tau[theta] == .(k_se) ~ "," ~ .(K_CORE) ~ "DiD core"),
    bquote(tau[theta] == .(s_se) ~ "," ~ .(S_CORE) ~ "DiD core")
  )
  gD <- ggplot2::ggplot(pd, ggplot2::aes(added, mean_posterior_sd,
                                         colour = setting, linetype = curve)) +
    ggplot2::geom_line(linewidth = 0.5) +
    # +/- 1.96 MC standard errors of the mean across replicates.
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = mean_posterior_sd - 1.96 * se_posterior_sd,
                   ymax = mean_posterior_sd + 1.96 * se_posterior_sd),
      size = 0.35, linewidth = 0.4, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c(low = "#0072B2", high = "#D55E00"),
                                 labels = pd_labs, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c("RCT + PP" = "solid", "DiD (reference)" = "dashed"),
      labels = c("RCT + PP" = "mixed: RCT + PP added",
                 "DiD (reference)" = "DiD added (reference)"),
      name = NULL) +
    ggplot2::labs(title = "B  Marginal value of incomplete designs (K, S)",
                  x = "Studies added to the DiD core",
                  y = expression("Posterior SD of " * mu[theta])) +
    .fig_theme()

  # --- I: efficiency over the trend-variability sweep (M) -----------------
  # Same sweep as panel C, on RMSE. mu_beta = 0 throughout, so no arm is
  # biased in mean and RMSE here is driven by the sampling VARIANCE of the
  # point estimate -- which is where the naive model's cost shows up, and
  # much more sharply than interval width shows it (naive/full RMSE reaches
  # ~5x on this sweep vs ~1.5x on width). Ignoring genuine between-study
  # trend variation propagates that variation straight into mu_theta.
  # metadid/naive RMSE ratio, paired per replicate where `paired` is supplied
  # (see paired_rmse_ratio()); otherwise falls back to the ratio of separately
  # aggregated RMSEs so a bare figure_panels.csv still renders locally.
  # One composition only (M: 15 DiD + 15 PP) -- the PP-heavy T sweep was not
  # resolvable at the replication counts in use, so it is not shown.
  f_ids <- scenario_ids("M")
  f_meta <- purrr::map_dfr(f_ids, function(id) {
    tibble::tibble(scenario_id = id,
                   x = SCENARIO_CONFIGS[[id]]$dgp$sigma_trend / 0.02)
  })
  pi_ <- if (!is.null(paired) && nrow(paired)) {
    paired |> dplyr::inner_join(f_meta, by = "scenario_id")
  } else {
    w <- te |>
      dplyr::filter(scenario_id %in% f_ids,
                    model_label %in% c("full", "naive")) |>
      dplyr::select(scenario_id, model_label, rmse) |>
      tidyr::pivot_wider(names_from = model_label, values_from = rmse)
    if (.has_arms(w)) {
      w |>
        dplyr::mutate(ratio = full / naive, lo = NA_real_, hi = NA_real_) |>
        dplyr::inner_join(f_meta, by = "scenario_id")
    } else {
      tibble::tibble(x = numeric(), ratio = numeric(), lo = numeric(),
                     hi = numeric())
    }
  }
  gI <- ggplot2::ggplot(pi_, ggplot2::aes(x, ratio)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "#6b6b6b") +
    ggplot2::geom_line(linewidth = 0.5, colour = "#0072B2") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi),
                             colour = "#0072B2", size = 0.35, linewidth = 0.4,
                             na.rm = TRUE) +
    ggplot2::labs(title = "D  Variable trends, zero mean: RMSE ratio (M)",
                  x = expression(tau[beta] * " (multiples of default), " *
                                   mu[beta] * " = 0"),
                  y = "RMSE of metadid / naive") +
    .fig_theme()

  # --- J: DiD studies needed to buy back the premium (R) ------------------
  # Ratio, not absolute width: total N grows along the sweep, so both models
  # tighten and only the RATIO isolates what the full model actually costs.
  # 1.0 means the insurance is free. One line per PP count, because the
  # premium is driven by how much PP evidence has to lean on the borrowed
  # trend -- not by the DiD count alone. Uses the paired per-replicate
  # contrast when supplied, so the pointranges are honest MC intervals.
  r_meta <- purrr::map_dfr(scenario_ids("R"), function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(scenario_id = id, x = d$n_did, n_pp = d$n_pp)
  })
  pj <- if (!is.null(paired) && nrow(paired)) {
    paired |> dplyr::inner_join(r_meta, by = "scenario_id")
  } else {
    w <- te |>
      dplyr::filter(scenario_id %in% r_meta$scenario_id,
                    model_label %in% c("full", "naive")) |>
      dplyr::select(scenario_id, model_label, mean_ci_width) |>
      tidyr::pivot_wider(names_from = model_label, values_from = mean_ci_width)
    if (.has_arms(w)) {
      w |>
        dplyr::mutate(ratio = full / naive, lo = NA_real_, hi = NA_real_) |>
        dplyr::inner_join(r_meta, by = "scenario_id")
    } else {
      tibble::tibble(x = numeric(), ratio = numeric(), lo = numeric(),
                     hi = numeric(), n_pp = numeric())
    }
  }
  pj <- dplyr::mutate(pj, n_pp = factor(n_pp, levels = sort(unique(n_pp))))
  gJ <- ggplot2::ggplot(pj, ggplot2::aes(x, ratio, colour = n_pp)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "#6b6b6b") +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi), size = 0.35,
                             linewidth = 0.4, na.rm = TRUE) +
    # ordered quantity -> sequential ramp, not categorical hues
    ggplot2::scale_colour_manual(values = c("#9ECAE1", "#4292C6", "#08519C"),
                                 name = "PP studies") +
    ggplot2::labs(title = "F  Pure null: interval premium (R)",
                  x = "DiD studies",
                  y = "90% CrI width, metadid / naive") +
    .fig_theme()

  # --- K/L: bias over the trend plane (Q) ---------------------------------
  # Each model is unbiased only on its own line: the full model along the
  # diagonal pp = did, the naive model along the horizontal pp = 0. The N and
  # O sweeps are a single vertical slice of this plane at did = -0.04.
  q_ids <- scenario_ids("Q")
  q_xy <- purrr::map_dfr(q_ids, function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(scenario_id = id,
                   x = d$did_trend / .fig_base(id),
                   y = d$pp_trend  / .fig_base(id))
  })
  pq <- te |>
    dplyr::filter(scenario_id %in% q_ids) |>
    dplyr::left_join(q_xy, by = "scenario_id")
  q_lim <- suppressWarnings(max(abs(pq$mean_bias), na.rm = TRUE))
  if (!is.finite(q_lim) || q_lim <= 0) q_lim <- 1  # empty/absent Q category

  .plane <- function(arm, ttl, slope) {
    ggplot2::ggplot(dplyr::filter(pq, model_label == arm),
                    ggplot2::aes(x, y, fill = mean_bias)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
      ggplot2::geom_abline(slope = slope, intercept = 0,
                           colour = "#111111", linewidth = 0.7) +
      .fig_bias_fill(q_lim) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = ttl,
                    x = "Mean time trend, DiD (normalised)",
                    y = "Mean time trend, pre-post (normalised)") +
      .fig_theme() +
      # .fig_theme() blanks legend titles, which suits categorical series but
      # leaves a continuous bias scale unlabelled -- restore it here.
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     legend.position = "right",
                     legend.title = ggplot2::element_text(size = 8))
  }
  gK <- .plane("full",  "G  Trend plane: metadid bias (Q)",  1)
  gL <- .plane("naive", "H  Trend plane: naive bias (Q)", 0)

  # --- M: which model to use, over the trend plane (Q) --------------------
  # K and L collapsed into a decision map. A cell is "neither usable" when
  # even the better of the two models is biased by more than `tol` -- there,
  # the PP studies carry no usable information about theta and the fallback
  # is a DiD-only analysis, which differences out its own trend and so is
  # unbiased anywhere on this plane.
  # Tolerable |bias| as a fraction of the true effect. This is a judgement
  # call, not a property of the data -- it is a plot-time argument so it can be
  # varied against an existing figure_panels.csv without rerunning anything.
  q_tol <- q_tol_frac * stats::median(abs(pq$true_value), na.rm = TRUE)
  pm <- pq |>
    dplyr::select(scenario_id, model_label, mean_bias, x, y) |>
    tidyr::pivot_wider(names_from = model_label, values_from = mean_bias)
  pm <- if (.has_arms(pm)) {
    # Four regions, not three. Where BOTH models are inside tolerance the
    # honest answer is that the choice does not matter -- collapsing that into
    # "metadid better" on a hair's-breadth difference overstates the case.
    dplyr::mutate(pm, region = dplyr::case_when(
      pmax(abs(full), abs(naive)) <= q_tol ~ "Either works",
      pmin(abs(full), abs(naive)) >  q_tol ~ "Neither usable: drop PP",
      abs(full) <= abs(naive)              ~ "metadid better",
      TRUE                                 ~ "naive better"
    ))
  } else {
    tibble::tibble(x = numeric(), y = numeric(), region = character())
  }
  gM <- ggplot2::ggplot(pm, ggplot2::aes(x, y, fill = region)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#111111",
                         linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, colour = "#111111", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = c(
      "Either works"            = "#009E73",
      "metadid better"          = "#0072B2",
      "naive better"            = "#D55E00",
      "Neither usable: drop PP" = "#4d4d4d"),
      breaks = c("Either works", "metadid better", "naive better",
                 "Neither usable: drop PP"), name = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "I  Which model to use (Q)",
      x = "Mean time trend, DiD (normalised)",
      y = "Mean time trend, pre-post (normalised)") +
    .fig_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 7,
                                                         colour = "#5c5c5c"))

  # --- U: what each incomplete design is worth, against rho ---------------
  # Information added PER STUDY, relative to one DiD study. Precision (not SD)
  # because information adds roughly linearly, which is what makes "per study"
  # meaningful. Thin dashed lines are the analytic per-study predictions
  #   PP / DiD = 2,   RCT / DiD = 2 (1 - rho)
  # so the vertical gap between a measured curve and its own prediction is the
  # price of incompleteness: PP must borrow the trend, the RCT must
  # accommodate baseline imbalance.
  u_ids     <- scenario_ids("U")
  u_meta    <- .u_meta(u_ids)
  u_anchors <- sort(unique(u_meta$anchor))
  pu_all    <- .u_measured(te, u_ids)
  u_th      <- .u_theory(u_ids, u_anchors)
  u_ceil    <- .u_ceiling(u_ids)

  # One panel per design, sharing construction. Both the Figure 1 panel and
  # the supplement are THEORY ONLY -- closed-form curves from equations (5.2)
  # and (6.5) of docs/design-information-derivation.docx. Set show_measured to
  # overlay the simulated points; nothing does at present.
  .exchange_panel <- function(designs, ttl, facet = FALSE,
                              show_measured = FALSE) {
    lab_anchor <- function(x) paste(x, "DiD anchor")
    th <- dplyr::filter(u_th, design %in% designs) |>
      dplyr::mutate(anchor = factor(lab_anchor(anchor)))
    ce <- dplyr::filter(u_ceil, design %in% designs)
    g <- ggplot2::ggplot(mapping = ggplot2::aes(rho, ratio)) +
      ggplot2::geom_hline(yintercept = 1, colour = "#dcdcdc") +
      # ceiling: what the design would be worth if the nuisance were KNOWN
      ggplot2::geom_line(data = ce, ggplot2::aes(linetype = "nuisance known"),
                         colour = "#6b6b6b", linewidth = 0.45) +
      # realised: the nuisance is estimated from the anchor
      ggplot2::geom_line(data = th,
                         ggplot2::aes(colour = anchor,
                                      linetype = "nuisance estimated"),
                         linewidth = 0.7)
    if (show_measured) {
      m <- dplyr::filter(pu_all, design %in% designs) |>
        dplyr::mutate(anchor = factor(lab_anchor(anchor)))
      g <- g +
        ggplot2::geom_point(data = m, ggplot2::aes(colour = anchor), size = 1.6,
                            na.rm = TRUE)
    }
    g <- g +
      ggplot2::scale_colour_manual(values = c("#9ECAE1", "#08519C"),
                                   name = NULL) +
      ggplot2::scale_linetype_manual(
        values = c("nuisance known" = "dotdash",
                   "nuisance estimated" = "solid"),
        breaks = c("nuisance known", "nuisance estimated"), name = NULL) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(title = ttl,
                    x = expression("Pre-post correlation " * rho),
                    y = "Information per study, relative to DiD") +
      .fig_theme()
    if (facet) g <- g + ggplot2::facet_wrap(~ design)
    g
  }
  gU  <- .exchange_panel("RCT (post-only)",
                         "E  Value of a post-only RCT, vs one DiD study")
  gU2 <- .exchange_panel(c("PP", "RCT (post-only)"),
                         "Value of each incomplete design, vs one DiD study",
                         facet = TRUE)

  list(A = gA,   # calibration (A/F)
       B = gB,   # trend-mean sweep, bias (J)
       C = gC,   # trend-variability sweep, coverage (M)
       D = gI,   # trend-variability sweep, RMSE (M)
       E = gD,   # composition sweeps, low vs high heterogeneity (K, S)
       F = gJ,   # pure-null interval premium (R)
       G = gK,   # trend plane, full model bias (Q)
       H = gL,   # trend plane, naive model bias (Q)
       I = gM,   # trend plane, which model to use (Q)
       U = gU,   # post-only RCT exchange rate vs rho (U)
       U2 = gU2) # both designs, faceted -- supplementary
}

illustration_panels <- function(draws) {
  truth <- draws$truth[1]
  # Labels are derived from the draws so the panel follows run_illustration()'s
  # n_studies rather than hard-coding 25/50.
  h <- if ("half" %in% names(draws)) draws$half[1] else 25L
  n <- if ("n_studies" %in% names(draws)) draws$n_studies[1] else 50L
  nm <- c(did_half_full  = sprintf("%d DiD alone", h),
          pp_half_naive  = sprintf("%d PP alone (no-trend)", h),
          mixed_full     = sprintf("metadid: %d DiD + %d PP", h, h),
          mixed_naive    = sprintf("naive: %d DiD + %d PP", h, h),
          did_all_full   = sprintf("all %d DiD (oracle)", n))
  lab <- nm
  cols <- setNames(c("#5b6770", "#D55E00", "#0072B2", "#D55E00", "#1a1a1a"),
                   unname(nm))
  d <- draws |> dplyr::mutate(model = unname(lab[model]))

  # Shared x range across BOTH panels, computed from every arm, so the split
  # and pooled panels can be compared directly rather than each self-scaling.
  xr <- range(d$draw, na.rm = TRUE)
  xr <- xr + c(-1, 1) * 0.04 * diff(xr)

  mk <- function(models, ttl, legend_rows = 1L) {
    dd <- d |> dplyr::filter(model %in% unname(lab[models])) |>
      dplyr::mutate(model = factor(model, levels = unname(lab[models])))
    ggplot2::ggplot(dd, ggplot2::aes(draw, colour = model, fill = model)) +
      ggplot2::geom_density(alpha = 0.25, adjust = 1.2, linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = truth, linetype = "dashed",
                          colour = "#1a1a1a") +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::coord_cartesian(xlim = xr) +
      # Three long series labels overflow a single row at panel width.
      ggplot2::guides(
        colour = ggplot2::guide_legend(nrow = legend_rows),
        fill   = ggplot2::guide_legend(nrow = legend_rows)) +
      ggplot2::labs(title = ttl,
                    x = expression("Population treatment effect " *
                                     mu[theta] * " (normalised)"),
                    y = NULL) +
      .fig_theme() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 7.5))
  }
  list(
    split  = mk(c("did_half_full", "pp_half_naive"),
                "One dataset, split and analysed separately"),
    # Three densities only: the two pooled analyses against the all-DiD
    # oracle. The DiD-half-alone arm carries panel A and is left out here so
    # the comparison that matters -- pooled vs oracle -- is not crowded.
    pooled = mk(c("mixed_full", "mixed_naive", "did_all_full"),
                "Pooling the designs", legend_rows = 2L)
  )
}

.retitle <- function(g, ttl) g + ggplot2::labs(title = ttl)

#' Figure 1: what the method does, the value of pooling, and its safety
#'
#' Column 1: the one-dataset illustration (split / pooled). Column 2: the
#' composition sweep at default and high effect heterogeneity. Column 3:
#' the pure-null interval premium and the zero-mean RMSE ratio.
#'
#' @param all_agg Aggregated results tibble.
#' @param illustration Draws tibble from run_illustration().
#' @return A patchwork object (6 panels, 2 rows x 3 columns).
plot_figure1 <- function(all_agg, illustration, paired = NULL,
                         q_tol_frac = Q_TOL_FRAC) {
  p  <- paper_panels(all_agg, paired, q_tol_frac)
  il <- illustration_panels(illustration)
  patchwork::wrap_plots(
    .retitle(il$split,  "A  Stratified inference"),
    .retitle(p$E, expression(bold("B  Marginal value of incomplete designs, low vs high ") * bold(tau[theta]))),
    .retitle(p$F, "C  Zero time trends"),
    .retitle(il$pooled, "D  Pooled inference"),
    .retitle(p$U, "E  Value of each design, vs one DiD study"),
    .retitle(p$D, "F  Mean zero time trends"),
    ncol = 3
  )
}

#' Supplementary figure: the exchange rate for both incomplete designs
#'
#' Figure 1 panel E shows only the post-only RCT. This keeps the PP comparison,
#' faceted, for the supplement.
#'
#' @param all_agg Aggregated results tibble.
#' @return A patchwork/ggplot object.
plot_figure_si <- function(all_agg) {
  paper_panels(all_agg)$U2
}

#' Figure 2: the trend plane and the model-choice map
#'
#' @param all_agg Aggregated results tibble.
#' @return A patchwork object (3 panels).
plot_figure2 <- function(all_agg, q_tol_frac = Q_TOL_FRAC) {
  p <- paper_panels(all_agg, q_tol_frac = q_tol_frac)
  patchwork::wrap_plots(
    .retitle(p$G, "A  metadid bias"),
    .retitle(p$H, "B  naive bias"),
    .retitle(p$I, "C  Which model to use"),
    ncol = 3
  ) +
    # A and B share an identical fill scale, so collecting guides merges them
    # into one legend and lets all three panels sit on a single row.
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(1.4, "lines"))
}

# ===========================================================================
# Figure rendering, decoupled from simulation
# ===========================================================================

#' Render Figure 1 and Figure 2 from exported CSVs
#'
#' The figures depend on the pipeline only through three CSVs, so plotting can
#' be iterated on without rerunning any simulation:
#'
#'   figure_panels.csv      per-scenario aggregates (incl. MC standard errors)
#'   illustration_draws.csv posterior draws for the column-1 densities
#'   paired_contrasts.csv   paired per-replicate ratios for panels C and F
#'
#' Download those from a `paper-figure` run artifact into `dir`, then:
#'
#'   for (f in c("R/scenarios.R", "R/assess.R", "R/plots.R")) source(f)
#'   render_paper_figures("output")
#'
#' figure.yaml calls this same function, so CI and local renders cannot drift.
#'
#' @param dir Directory holding the CSVs.
#' @param out Directory to write figures to (defaults to `dir`).
#' @param save Write PNG/PDF files; FALSE just returns the plot objects.
#' @return Invisibly, a list of the two patchwork objects.
render_paper_figures <- function(dir = "output", out = dir, save = TRUE,
                                 q_tol_frac = Q_TOL_FRAC) {
  rd <- function(f, required = TRUE) {
    path <- file.path(dir, f)
    if (!file.exists(path)) {
      if (required) stop("missing required input: ", path, call. = FALSE)
      return(NULL)
    }
    readr::read_csv(path, show_col_types = FALSE)
  }
  agg    <- rd("figure_panels.csv")
  illus  <- rd("illustration_draws.csv")
  paired <- rd("paired_contrasts.csv", required = FALSE)

  figs <- list(
    figure1 = list(plot = plot_figure1(agg, illus, paired, q_tol_frac),
                   w = 13, h = 9),
    figure2 = list(plot = plot_figure2(agg, q_tol_frac), w = 13, h = 5),
    figure_si = list(plot = plot_figure_si(agg), w = 9, h = 4.6),
    # theory vs simulation, one per incomplete design
    figure_val_pp  = list(plot = plot_design_validation(agg, "PP"),
                          w = 10, h = 8),
    figure_val_rct = list(plot = plot_design_validation(agg, "RCT (post-only)"),
                          w = 10, h = 8)
  )
  if (save) {
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    for (nm in names(figs)) {
      f <- figs[[nm]]
      ggplot2::ggsave(file.path(out, paste0(nm, ".png")), f$plot,
                      width = f$w, height = f$h, dpi = 300, limitsize = FALSE)
      ggplot2::ggsave(file.path(out, paste0(nm, ".pdf")), f$plot,
                      width = f$w, height = f$h, limitsize = FALSE)
    }
  }
  invisible(lapply(figs, `[[`, "plot"))
}

#' How the panel-I regions respond to the bias tolerance
#'
#' The "which model to use" map depends on a judgement call -- how much bias in
#' mu_theta is tolerable -- not on the data alone. This tabulates the region
#' counts across candidate tolerances so the choice can be made on evidence and
#' stated in the caption, rather than left implicit at whatever the default is.
#' Runs entirely from an existing figure_panels.csv; no simulation needed.
#'
#' @param all_agg Aggregated results (or the contents of figure_panels.csv).
#' @param fracs Candidate values of q_tol_frac.
#' @return Tibble: one row per tolerance, one column per region, plus the
#'   absolute bias threshold that tolerance implies.
q_tolerance_scan <- function(all_agg, fracs = c(0.05, 0.10, 0.15, 0.20,
                                                0.25, 0.33, 0.50)) {
  purrr::map_dfr(fracs, function(f) {
    d <- paper_panels(all_agg, q_tol_frac = f)$I$data
    if (!nrow(d)) return(tibble::tibble(q_tol_frac = f))
    truth <- stats::median(abs(all_agg$true_value), na.rm = TRUE)
    tibble::tibble(q_tol_frac = f, abs_threshold = f * truth) |>
      dplyr::bind_cols(
        as.list(table(factor(d$region,
          levels = c("Either works", "metadid better", "naive better",
                     "Neither usable: drop PP"))))
      )
  })
}

#' Theory-vs-empirical validation figure for one incomplete design
#'
#' Four panels:
#'   1-2  exchange rate against rho, one per anchor: the closed-form curve
#'        (equation 6.5) with the simulated points on top
#'   3    the absolute posterior SDs the ratio is built from, so any
#'        disagreement can be traced to a particular arm rather than guessed at
#'   4    predicted against observed, with a 1:1 line
#'
#' Panels 1-2 and 4 answer "do they agree"; panel 3 answers "and if not, where
#' does it come from".
#'
#' @param all_agg Aggregated results tibble (or figure_panels.csv).
#' @param design  "PP" or "RCT (post-only)".
#' @return A patchwork object (4 panels, 2 x 2).
plot_design_validation <- function(all_agg, design = "PP") {
  te <- dplyr::filter(all_agg, parameter == "treatment_effect_mean")
  ids <- scenario_ids("U")
  meta <- .u_meta(ids)
  anchors <- sort(unique(meta$anchor))
  rhos <- sort(unique(meta$rho))
  arm_key <- if (design == "PP") "pp" else "rct"

  meas <- dplyr::filter(.u_measured(te, ids), design == !!design)
  th   <- dplyr::filter(.u_theory(ids, anchors), design == !!design)
  # theory evaluated AT the simulated rho values, for panel 4
  th_at <- dplyr::filter(.u_theory(ids, anchors, rho = rhos),
                         design == !!design) |>
    dplyr::select(rho, anchor, predicted = ratio)

  by_anchor <- function(a) {
    ggplot2::ggplot(mapping = ggplot2::aes(rho, ratio)) +
      ggplot2::geom_hline(yintercept = 1, colour = "#dcdcdc") +
      ggplot2::geom_line(data = dplyr::filter(th, anchor == a),
                         colour = "#08519C", linewidth = 0.6) +
      ggplot2::geom_point(data = dplyr::filter(meas, anchor == a),
                          colour = "#111111", size = 2, na.rm = TRUE) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(title = sprintf("%d DiD anchor", a),
                    x = expression(rho),
                    y = "Information per study, rel. DiD") +
      .fig_theme()
  }

  # panel 3: the absolute quantities the ratio is derived from
  abs_dat <- te |>
    dplyr::filter(scenario_id %in% ids) |>
    dplyr::inner_join(meta, by = "scenario_id") |>
    dplyr::filter(arm %in% c("core", "did", arm_key)) |>
    dplyr::mutate(arm = factor(arm, levels = c("core", "did", arm_key),
                               labels = c("core (DiD only)", "+ DiD",
                                          paste("+", design))),
                  anchor = factor(paste(anchor, "anchor")))
  g3 <- ggplot2::ggplot(abs_dat, ggplot2::aes(rho, mean_posterior_sd,
                                              colour = arm, linetype = anchor)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = c("#6b6b6b", "#08519C", "#CC79A7"),
                                 name = NULL) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(title = "Absolute inputs",
                  x = expression(rho),
                  y = expression("Posterior SD of " * mu[theta])) +
    .fig_theme()

  pvo <- meas |>
    dplyr::inner_join(th_at, by = c("rho", "anchor")) |>
    dplyr::mutate(anchor = factor(paste(anchor, "anchor")))
  lim <- suppressWarnings(range(c(pvo$ratio, pvo$predicted, 1), na.rm = TRUE))
  if (!all(is.finite(lim)) || diff(lim) == 0) lim <- c(0, 1.2)
  g4 <- ggplot2::ggplot(pvo, ggplot2::aes(predicted, ratio, colour = anchor)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#6b6b6b",
                         linetype = "dashed") +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c("#9ECAE1", "#08519C"), name = NULL) +
    ggplot2::coord_fixed(xlim = lim, ylim = lim) +
    ggplot2::labs(title = "Predicted vs observed",
                  x = "Predicted", y = "Observed") +
    .fig_theme()

  patchwork::wrap_plots(by_anchor(anchors[1]), by_anchor(anchors[2]), g3, g4,
                        ncol = 2) +
    patchwork::plot_annotation(
      title = paste0("Theory vs simulation: ", design),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11)))
}
