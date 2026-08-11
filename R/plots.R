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
#' @param all_agg Aggregated results tibble (scenario x model x parameter).
#' @return Named list of ggplot objects, one per panel letter.
#'
#' Split out from plot_paper_figure() so a single panel can be rebuilt and
#' previewed locally from a downloaded figure_panels.csv, without rerunning
#' the pipeline:  paper_panels(readr::read_csv("figure_panels.csv"))$I
#' Panels whose category is absent from `all_agg` come back empty rather than
#' erroring, so a partial CSV still lets you inspect the panels it does cover.
paper_panels <- function(all_agg, paired = NULL) {
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

  # --- D: information from incomplete designs (K) -------------------------
  k_ids <- scenario_ids("K")
  k_info <- purrr::map_dfr(k_ids, function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(
      scenario_id = id,
      added = (d$n_did - K_CORE) + d$n_rct + d$n_pp,
      curve = dplyr::case_when(
        d$n_rct > 0 & d$n_pp > 0 ~ "RCT + PP",
        d$n_rct > 0              ~ "RCT only",
        d$n_pp  > 0              ~ "PP only",
        d$n_did > K_CORE         ~ "DiD (reference)",
        TRUE                     ~ "core"
      )
    )
  })
  # Only the mixed RCT+PP sweep against the DiD reference, matching the S
  # panel so the two are read side by side. The RCT-only and PP-only arms are
  # still simulated and remain available for the SI.
  k_info <- dplyr::filter(k_info,
                          curve %in% c("RCT + PP", "DiD (reference)", "core"))
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
    ggplot2::labs(title = "E  Adding incomplete designs to a small DiD core (K)",
                  x = sprintf("Studies added to a core of %d DiD", K_CORE),
                  y = expression("Posterior SD of " * mu[theta])) +
    ggplot2::expand_limits(y = 0) +
    .fig_theme()

  # --- I: efficiency over the trend-variability sweep (M) -----------------
  # Same sweep as panel C, on RMSE. mu_beta = 0 throughout, so no arm is
  # biased in mean and RMSE here is driven by the sampling VARIANCE of the
  # point estimate -- which is where the naive model's cost shows up, and
  # much more sharply than interval width shows it (naive/full RMSE reaches
  # ~5x on this sweep vs ~1.5x on width). Ignoring genuine between-study
  # trend variation propagates that variation straight into mu_theta.
  # Full/naive RMSE ratio, paired per replicate where `paired` is supplied
  # (see paired_rmse_ratio()); otherwise falls back to the ratio of separately
  # aggregated RMSEs so a bare figure_panels.csv still renders locally.
  # Two lines: M is 15 DiD + 15 PP, T is 5 DiD + 25 PP. The more the evidence
  # base leans on PP studies, the earlier ignoring trend variation costs you.
  f_ids <- c(scenario_ids("M"), scenario_ids("T"))
  f_meta <- purrr::map_dfr(f_ids, function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(scenario_id = id,
                   x = d$sigma_trend / 0.02,
                   composition = sprintf("%d DiD + %d PP", d$n_did, d$n_pp))
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
                     hi = numeric(), composition = character())
    }
  }
  gI <- ggplot2::ggplot(pi_, ggplot2::aes(x, ratio, colour = composition)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "#6b6b6b") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = composition),
                         alpha = 0.15, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = c("#0072B2", "#CC79A7"),
                                 name = NULL, aesthetics = c("colour", "fill")) +
    ggplot2::labs(title = "D  Variable trends, zero mean: RMSE ratio (M, T)",
                  x = expression(tau[beta] * " (multiples of default), " *
                                   mu[beta] * " = 0"),
                  y = "RMSE of metadid / naive") +
    .fig_theme()

  # --- J: DiD studies needed to buy back the premium (R) ------------------
  # Ratio, not absolute width: total N grows along the sweep, so both models
  # tighten and only the RATIO isolates what the full model actually costs.
  # 1.0 means the insurance is free. One line per PP count, because the
  # premium is driven by how much PP evidence has to lean on the borrowed
  # trend -- not by the DiD count alone.
  r_meta <- purrr::map_dfr(scenario_ids("R"), function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(scenario_id = id, x = d$n_did, n_pp = d$n_pp)
  })
  pj <- te |>
    dplyr::filter(scenario_id %in% r_meta$scenario_id,
                  model_label %in% c("full", "naive")) |>
    dplyr::select(scenario_id, model_label, mean_ci_width) |>
    tidyr::pivot_wider(names_from = model_label, values_from = mean_ci_width)
  pj <- if (.has_arms(pj)) {
    pj |>
      dplyr::inner_join(r_meta, by = "scenario_id") |>
      dplyr::mutate(ratio = full / naive,
                    n_pp = factor(n_pp, levels = sort(unique(n_pp))))
  } else {
    tibble::tibble(x = numeric(), ratio = numeric(), n_pp = factor())
  }
  gJ <- ggplot2::ggplot(pj, ggplot2::aes(x, ratio, colour = n_pp)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "#6b6b6b") +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.6) +
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
  q_tol_frac <- 0.10   # tolerable |bias| as a fraction of the true effect
  q_tol <- q_tol_frac * stats::median(abs(pq$true_value), na.rm = TRUE)
  pm <- pq |>
    dplyr::select(scenario_id, model_label, mean_bias, x, y) |>
    tidyr::pivot_wider(names_from = model_label, values_from = mean_bias)
  pm <- if (.has_arms(pm)) {
    dplyr::mutate(pm, region = dplyr::case_when(
      pmin(abs(full), abs(naive)) > q_tol ~ "Neither usable: drop PP",
      abs(full) <= abs(naive)             ~ "metadid better",
      TRUE                                ~ "naive better"
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
      "metadid better"          = "#0072B2",
      "naive better"            = "#D55E00",
      "Neither usable: drop PP" = "#4d4d4d"), name = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "I  Which model to use (Q)",
      subtitle = sprintf("'neither' = best available |bias| > %.0f%% of the true effect",
                         100 * q_tol_frac),
      x = "Mean time trend, DiD (normalised)",
      y = "Mean time trend, pre-post (normalised)") +
    .fig_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 7,
                                                         colour = "#5c5c5c"))

  # --- S: composition sweep under high effect heterogeneity (S) -----------
  # Same construction as the K panel, on the S grid (sigma_effect = 0.10).
  # Only two curves exist here: the mixed RCT+PP sweep and the DiD reference.
  s_ids <- scenario_ids("S")
  s_info <- purrr::map_dfr(s_ids, function(id) {
    d <- SCENARIO_CONFIGS[[id]]$dgp
    tibble::tibble(
      scenario_id = id,
      added = (d$n_did - S_CORE) + d$n_rct + d$n_pp,
      curve = dplyr::case_when(
        d$n_rct > 0 & d$n_pp > 0 ~ "RCT + PP",
        d$n_did > S_CORE         ~ "DiD (reference)",
        TRUE                     ~ "core"
      )
    )
  })
  s_curves <- setdiff(unique(s_info$curve), "core")
  s_core <- s_info |> dplyr::filter(curve == "core")
  s_all <- dplyr::bind_rows(
    s_info |> dplyr::filter(curve != "core"),
    tidyr::crossing(s_core |> dplyr::select(-curve), curve = s_curves)
  )
  ps <- te |>
    dplyr::filter(scenario_id %in% s_ids) |>
    dplyr::inner_join(s_all, by = "scenario_id",
                      relationship = "many-to-many")
  gS <- ggplot2::ggplot(ps, ggplot2::aes(added, mean_posterior_sd,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = c(
      "RCT + PP" = "#0072B2", "DiD (reference)" = "#6b6b6b")) +
    ggplot2::scale_linetype_manual(values = c(
      "RCT + PP" = "solid", "DiD (reference)" = "dashed")) +
    ggplot2::labs(
      title = "Incomplete designs reach parity under high heterogeneity (S)",
      x = sprintf("Studies added to a core of %d DiD", S_CORE),
      y = expression("Posterior SD of " * mu[theta])) +
    ggplot2::expand_limits(y = 0) +
    .fig_theme()

  list(A = gA,   # calibration (A/F)
       B = gB,   # trend-mean sweep, bias (J)
       C = gC,   # trend-variability sweep, coverage (M)
       D = gI,   # trend-variability sweep, RMSE (M)
       E = gD,   # information from incomplete designs (K)
       F = gJ,   # pure-null interval premium (R)
       G = gK,   # trend plane, full model bias (Q)
       H = gL,   # trend plane, naive model bias (Q)
       I = gM,   # trend plane, which model to use (Q)
       S = gS)   # composition sweep, high effect heterogeneity (S)
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

  mk <- function(models, ttl) {
    dd <- d |> dplyr::filter(model %in% unname(lab[models])) |>
      dplyr::mutate(model = factor(model, levels = unname(lab[models])))
    ggplot2::ggplot(dd, ggplot2::aes(draw, colour = model, fill = model)) +
      ggplot2::geom_density(alpha = 0.25, adjust = 1.2, linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = truth, linetype = "dashed",
                          colour = "#1a1a1a") +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(title = ttl,
                    x = expression("Population treatment effect " *
                                     mu[theta] * " (normalised)"),
                    y = NULL) +
      .fig_theme() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }
  list(
    split  = mk(c("did_half_full", "pp_half_naive"),
                "One dataset, split and analysed separately"),
    # Three densities only: the two pooled analyses against the all-DiD
    # oracle. The DiD-half-alone arm carries panel A and is left out here so
    # the comparison that matters -- pooled vs oracle -- is not crowded.
    pooled = mk(c("mixed_full", "mixed_naive", "did_all_full"),
                "Pooling the designs")
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
plot_figure1 <- function(all_agg, illustration, paired = NULL) {
  p  <- paper_panels(all_agg, paired)
  il <- illustration_panels(illustration)
  patchwork::wrap_plots(
    .retitle(il$split,  "A  Stratified inference"),
    .retitle(p$E, expression("B  Marginal study information for low " * tau[theta])),
    .retitle(p$F, "C  Zero time trends"),
    .retitle(il$pooled, "D  Pooled inference"),
    .retitle(p$S, expression("E  Marginal study information for high " * tau[theta])),
    .retitle(p$D, "F  Mean zero time trends"),
    ncol = 3
  )
}

#' Figure 2: the trend plane and the model-choice map
#'
#' @param all_agg Aggregated results tibble.
#' @return A patchwork object (3 panels).
plot_figure2 <- function(all_agg) {
  p <- paper_panels(all_agg)
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
