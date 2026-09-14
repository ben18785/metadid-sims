# Performance over time, from archive.csv.
#
# Kept out of plots.R on purpose: three of these four functions are data
# functions, not plotting, and the unit suite skips plots.R to stay fast and
# dependency-light. Living here means they are covered by tests.

# ===========================================================================
# Performance over time (archive.csv)
# ===========================================================================
#
# archive.csv accumulates one row per (scenario, model, parameter) per run,
# stamped with run_date and metadid_install_datetime. It has been collecting
# since 2026-05 and nothing plots it.
#
# DESIGN CONSTRAINT: none of this may fail the pipeline. The archive is a
# nice-to-have diagnostic restored from gh-pages by a curl that is allowed to
# miss; it can be absent, truncated, malformed, or too short to trend. Every
# function here returns NULL or an empty tibble rather than raising, and the
# report chunks degrade to a note. A broken history panel must never cost a
# run that took hours to fit.

# Runs with very few replications are debug/smoke runs (one archived run used
# n_reps 3 and 6). Their metrics are not comparable and would read as spikes.
.ARCHIVE_MIN_REPS <- 10L

#' Load archive.csv defensively
#'
#' @return A tibble, or NULL if the archive is missing/unusable for any reason.
#'   Never raises.
load_archive <- function(path = "output/archive.csv") {
  tryCatch({
    if (!file.exists(path)) return(NULL)
    a <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    need <- c("scenario_id", "model_label", "parameter", "run_date", "n_reps",
              "empirical_coverage", "mean_bias", "rmse")
    if (!all(need %in% names(a))) return(NULL)
    a <- a[!is.na(a$run_date) & !is.na(a$n_reps), , drop = FALSE]
    # Drop smoke runs: judged per RUN, not per row, so a run is kept or
    # dropped whole and categories with legitimately smaller counts survive.
    keep <- a |>
      dplyr::group_by(.data$run_date) |>
      dplyr::summarise(ok = max(.data$n_reps, na.rm = TRUE) >= .ARCHIVE_MIN_REPS,
                       .groups = "drop") |>
      dplyr::filter(.data$ok) |>
      dplyr::pull("run_date")
    a <- a[a$run_date %in% keep, , drop = FALSE]
    if (nrow(a) == 0 || dplyr::n_distinct(a$run_date) < 2) return(NULL)
    a
  }, error = function(e) NULL, warning = function(w) NULL)
}

#' Restrict the archive to series with a comparable history
#'
#' The scenario set grew from 64 to 160 over the archived period, so any
#' aggregate over "all scenarios" moves with composition rather than
#' performance. Trends are therefore computed only on series present in at
#' least `min_runs` of the runs.
archive_common_series <- function(archive, min_runs = NULL) {
  if (is.null(archive)) return(NULL)
  tryCatch({
    n_runs <- dplyr::n_distinct(archive$run_date)
    if (is.null(min_runs)) min_runs <- max(2L, ceiling(n_runs * 0.75))
    counts <- archive |>
      dplyr::count(.data$scenario_id, .data$model_label, .data$parameter,
                   name = "n_runs_seen")
    keep <- counts[counts$n_runs_seen >= min_runs, , drop = FALSE]
    if (nrow(keep) == 0) return(NULL)
    dplyr::inner_join(
      archive, keep[, c("scenario_id", "model_label", "parameter")],
      by = c("scenario_id", "model_label", "parameter")
    )
  }, error = function(e) NULL)
}

#' Flag series whose latest value moved beyond ordinary run-to-run variation
#'
#' Compares the most recent run against the median of the preceding ones, in
#' units of that series' own historical spread (MAD). This is a CHANGE
#' DETECTOR, not a trend estimator: with a handful of irregularly spaced runs
#' it can catch a step change from a code change, and should not be read as
#' evidence about gradual drift.
#'
#' Most series are BIT-IDENTICAL between runs, because targets caches fits and
#' only rebuilds what changed. Their MAD is then zero, or floating-point noise
#' a few times 1e-17, and dividing by it turns a change in the fourth decimal
#' into a z-score of several hundred. A statistical criterion alone is
#' therefore meaningless here: the first version of this flagged 70 "rmse
#' regressions", essentially all of them cache artefacts.
#'
#' So a series must clear BOTH bars: unusual relative to its own history
#' (|z| >= `z`) AND materially changed (`min_rel_change`, default 5% of the
#' baseline). Materiality is what does the real work.
#'
#' @return A tibble of flagged series, possibly empty. Never raises.
#' `scale_metric` sets what "materially changed" is measured against. For rmse
#' and coverage the metric's own baseline is the right denominator. For
#' mean_bias it is not: bias sits near zero, so a move from -2e-04 to +1e-04 is
#' a "157% change" and means nothing. Pass scale_metric = "rmse" there, and a
#' bias shift is judged against the estimator's own error scale, which is the
#' quantity that decides whether it matters.
detect_regressions <- function(archive, metric = "rmse", z = 3,
                               min_rel_change = 0.05,
                               scale_metric = metric) {
  if (is.null(archive)) return(tibble::tibble())
  tryCatch({
    a <- archive[!is.na(archive[[metric]]), , drop = FALSE]
    if (nrow(a) == 0) return(tibble::tibble())
    latest <- max(a$run_date)
    a |>
      dplyr::group_by(.data$scenario_id, .data$model_label, .data$parameter) |>
      dplyr::filter(dplyr::n() >= 4, any(.data$run_date == latest)) |>
      dplyr::summarise(
        current  = .data[[metric]][.data$run_date == latest][1],
        baseline = stats::median(.data[[metric]][.data$run_date != latest]),
        spread   = stats::mad(.data[[metric]][.data$run_date != latest]),
        scale    = stats::median(.data[[scale_metric]][.data$run_date != latest],
                                 na.rm = TRUE),
        .groups  = "drop"
      ) |>
      dplyr::mutate(
        # Floor the spread at a fraction of the baseline's own magnitude. A MAD
        # below that is a cached, unchanged series, not a precisely-measured
        # one, and must not license a huge z-score.
        spread = ifelse(.data$spread <= abs(.data$baseline) * 1e-6,
                        NA_real_, .data$spread),
        change     = .data$current - .data$baseline,
        rel_change = ifelse(is.na(.data$scale) | .data$scale == 0, NA_real_,
                            .data$change / abs(.data$scale)),
        z_score    = ifelse(is.na(.data$spread), NA_real_,
                            .data$change / .data$spread)
      ) |>
      dplyr::filter(
        !is.na(.data$z_score), abs(.data$z_score) >= z,
        !is.na(.data$rel_change), abs(.data$rel_change) >= min_rel_change
      ) |>
      dplyr::arrange(dplyr::desc(abs(.data$rel_change)))
  }, error = function(e) tibble::tibble())
}

#' Trend plot of a metric over runs, faceted by category
#'
#' @return A ggplot, or NULL if there is nothing plottable. Never raises.
plot_performance_trend <- function(archive, metric = "rmse",
                                   parameter = "treatment_effect_mean") {
  if (is.null(archive)) return(NULL)
  tryCatch({
    d <- archive[archive$parameter == parameter, , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$category <- substr(d$scenario_id, 1L, 1L)
    d$series <- paste(d$scenario_id, d$model_label)
    # Version boundaries: where the installed metadid changed between runs.
    vlines <- NULL
    if ("metadid_install_datetime" %in% names(d)) {
      vr <- d |>
        dplyr::distinct(.data$run_date, .data$metadid_install_datetime) |>
        dplyr::arrange(.data$run_date)
      if (nrow(vr) > 1) {
        changed <- vr$metadid_install_datetime[-1] != vr$metadid_install_datetime[-nrow(vr)]
        vlines <- vr$run_date[-1][which(changed)]
      }
    }
    g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$run_date, y = .data[[metric]],
                                         group = .data$series)) +
      ggplot2::geom_line(alpha = 0.35, linewidth = 0.3) +
      ggplot2::geom_point(size = 0.6, alpha = 0.5) +
      ggplot2::facet_wrap(~ category, scales = "free_y") +
      ggplot2::labs(
        x = NULL, y = metric,
        title = paste0(metric, " over runs (", parameter, ")"),
        subtitle = paste("Series present in most runs only; dashed rules mark",
                         "a change in the installed metadid")
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (length(vlines)) {
      g <- g + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed",
                                   colour = "grey40", linewidth = 0.3)
    }
    g
  }, error = function(e) NULL)
}
