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
