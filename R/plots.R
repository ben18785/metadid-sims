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
  normal      = "#0072B2", robust = "#009E73"
)

#' Multi-panel paper figure from aggregated results
#'
#' @param all_agg Aggregated results tibble (scenario x model x parameter).
#' @return A patchwork object (8 panels, 4 x 2).
plot_paper_figure <- function(all_agg) {
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
  band <- 1.96 * sqrt(0.9 * 0.1 / max(pa$n_reps, na.rm = TRUE))
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

  # --- D: information from incomplete designs (K) -------------------------
  k_ids <- scenario_ids("K")
  k_info <- purrr::map_dfr(k_ids, function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(
      scenario_id = id,
      added = (d$n_did - 10L) + d$n_rct + d$n_pp,
      curve = dplyr::case_when(
        d$n_rct > 0 & d$n_pp > 0 ~ "RCT + PP",
        d$n_rct > 0              ~ "RCT only",
        d$n_pp  > 0              ~ "PP only",
        d$n_did > 10L            ~ "DiD (reference)",
        TRUE                     ~ "core"
      )
    )
  })
  curves <- setdiff(unique(k_info$curve), "core")
  k_core <- k_info |> dplyr::filter(curve == "core")
  k_all <- dplyr::bind_rows(
    k_info |> dplyr::filter(curve != "core"),
    tidyr::crossing(k_core |> dplyr::select(-curve), curve = curves)
  )
  pd <- te |>
    dplyr::filter(scenario_id %in% k_ids) |>
    dplyr::inner_join(k_all, by = "scenario_id",
                      relationship = "many-to-many")
  gD <- ggplot2::ggplot(pd, ggplot2::aes(added, mean_posterior_sd,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = c(
      "RCT + PP" = "#0072B2", "RCT only" = "#009E73",
      "PP only" = "#CC79A7", "DiD (reference)" = "#6b6b6b")) +
    ggplot2::scale_linetype_manual(values = c(
      "RCT + PP" = "solid", "RCT only" = "solid",
      "PP only" = "solid", "DiD (reference)" = "dashed")) +
    ggplot2::labs(title = "D  Information from incomplete designs (K)",
                  x = "Studies added to a core of 10 DiD",
                  y = expression("Posterior SD of " * mu[theta])) +
    ggplot2::expand_limits(y = 0) +
    .fig_theme()

  # --- E: outlier contamination (L) ---------------------------------------
  pe <- sweep_data("L", function(d, id) 100 * d$n_outlier / d$n_did)
  gE <- ggplot2::ggplot(pe, ggplot2::aes(x, rmse, colour = model_label)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "E  Outlying studies: RMSE (L)",
                  x = "Outlying studies (%)",
                  y = expression("RMSE of " * mu[theta])) +
    ggplot2::expand_limits(y = 0) +
    .fig_theme()

  # --- F: exchangeability gap sweep (N) -----------------------------------
  pf <- sweep_data("N", function(d, id) (d$did_trend - d$pp_trend) / .fig_base(id))
  gF <- ggplot2::ggplot(pf, ggplot2::aes(x, mean_bias, colour = model_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                         linetype = "dotted", formula = y ~ x) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "F  Exchangeability violations (N)",
                  x = expression("PP vs DiD trend difference " * Delta *
                                   " (normalised)"),
                  y = expression("Bias in " * mu[theta])) +
    .fig_theme()

  # --- G: crossover sweep (O) ---------------------------------------------
  pg <- sweep_data("O", function(d, id) d$pp_trend / .fig_base(id))
  x_did <- SCENARIO_CONFIGS[["O1"]]$dgp$did_trend / .fig_base("O1")
  gG <- ggplot2::ggplot(pg, ggplot2::aes(x, mean_bias, colour = model_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = x_did, linetype = "dotdash",
                        colour = "#6b6b6b") +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                         linetype = "dotted", formula = y ~ x) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "G  Trend crossover: when borrowing hurts (O)",
                  x = expression("True PP trend " * mu[beta]^{PP} *
                                   " (normalised); vline = DiD trend"),
                  y = expression("Bias in " * mu[theta])) +
    .fig_theme()

  # --- H: RCT baseline-imbalance sweep (P) --------------------------------
  ph <- sweep_data("P", function(d, id) d$rct_gamma_mean / .fig_base(id))
  gH <- ggplot2::ggplot(ph, ggplot2::aes(x, mean_bias, colour = model_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                         linetype = "dotted", formula = y ~ x) +
    ggplot2::scale_colour_manual(values = .fig_model_cols) +
    ggplot2::labs(title = "H  RCT baseline imbalance (P)",
                  x = expression("RCT imbalance " * mu[gamma] * " (normalised)"),
                  y = expression("Bias in " * mu[theta])) +
    .fig_theme()

  patchwork::wrap_plots(gA, gB, gC, gD, gE, gF, gG, gH, ncol = 2)
}
