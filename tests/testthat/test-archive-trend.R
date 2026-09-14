# The performance-over-time panel reads archive.csv, which is restored from
# gh-pages by a curl that is ALLOWED to miss. It can be absent, truncated,
# malformed, or too short to trend.
#
# The requirement is that none of that fails the pipeline. A run costs hours
# of MCMC; a broken history panel must never take one down. These tests feed
# the loader every broken input we can construct and assert it degrades to
# NULL or an empty tibble rather than raising.

suppressMessages(library(ggplot2))

tmp_csv <- function(txt) {
  p <- tempfile(fileext = ".csv"); writeLines(txt, p); p
}

good_archive <- function(n_runs = 5, n_series = 4) {
  runs <- sprintf("2026-0%d-01T00:00:00Z", seq_len(n_runs))
  tidyr::crossing(run_date = runs, scenario_id = paste0("A", seq_len(n_series))) |>
    dplyr::mutate(
      model_label = "default", parameter = "treatment_effect_mean",
      n_reps = 25L, empirical_coverage = 0.94,
      mean_bias = 0.001, rmse = 0.02,
      metadid_install_datetime = "2026-01-01T00:00:00Z"
    )
}

# --- the loader must never raise ------------------------------------------

test_that("a missing archive yields NULL, not an error", {
  expect_null(load_archive(file.path(tempdir(), "definitely-not-here.csv")))
})

test_that("an empty or header-only archive yields NULL", {
  expect_null(load_archive(tmp_csv("")))
  expect_null(load_archive(tmp_csv("scenario_id,run_date,n_reps")))
})

test_that("an archive missing required columns yields NULL", {
  expect_null(load_archive(tmp_csv(c("a,b,c", "1,2,3"))))
})

test_that("malformed, ragged and binary content yields NULL rather than raising", {
  expect_null(load_archive(tmp_csv(c("scenario_id,run_date", "A1"))))
  expect_null(load_archive(tmp_csv(c("<<<not a csv>>>", "\001\002\003"))))
})

test_that("a single-run archive yields NULL — nothing to trend", {
  p <- tempfile(fileext = ".csv")
  readr::write_csv(good_archive(n_runs = 1), p)
  expect_null(load_archive(p))
})

test_that("smoke runs are dropped whole, by run rather than by row", {
  # One archived run used n_reps 3 and 6. Its metrics are not comparable and
  # would read as a spike in every panel.
  a <- dplyr::bind_rows(
    good_archive(n_runs = 3),
    good_archive(n_runs = 1) |> dplyr::mutate(run_date = "2026-09-09T00:00:00Z",
                                              n_reps = 3L)
  )
  p <- tempfile(fileext = ".csv"); readr::write_csv(a, p)
  out <- load_archive(p)
  expect_false("2026-09-09T00:00:00Z" %in% out$run_date)
  expect_equal(dplyr::n_distinct(out$run_date), 3L)
})

# --- the consumers must never raise ---------------------------------------

test_that("every consumer tolerates NULL input", {
  expect_null(archive_common_series(NULL))
  expect_null(plot_performance_trend(NULL))
  expect_equal(nrow(detect_regressions(NULL)), 0L)
})

test_that("consumers tolerate an archive with no matching rows", {
  a <- good_archive()
  expect_null(plot_performance_trend(a, parameter = "no_such_parameter"))
  expect_equal(nrow(detect_regressions(a, metric = "rmse", z = 99)), 0L)
})

test_that("detect_regressions survives all-NA and zero-variance series", {
  a <- good_archive()                      # rmse identical in every run
  expect_equal(nrow(detect_regressions(a)), 0L)   # zero MAD must not divide by zero
  a$rmse <- NA_real_
  expect_equal(nrow(detect_regressions(a)), 0L)
})

# --- and it must actually work ---------------------------------------------

test_that("a genuine step change is detected", {
  a <- good_archive(n_runs = 6, n_series = 2)
  a$rmse <- a$rmse + stats::runif(nrow(a), -0.0005, 0.0005)   # ordinary jitter
  latest <- max(a$run_date)
  a$rmse[a$run_date == latest & a$scenario_id == "A1"] <- 0.20  # a real jump
  hits <- detect_regressions(a)
  expect_true("A1" %in% hits$scenario_id)
  expect_false("A2" %in% hits$scenario_id)
  expect_true(all(abs(hits$z_score) >= 3))
})

test_that("common-series filtering excludes scenarios added part-way", {
  a <- dplyr::bind_rows(
    good_archive(n_runs = 4, n_series = 2),
    good_archive(n_runs = 1, n_series = 1) |>
      dplyr::mutate(run_date = "2026-05-01T00:00:00Z", scenario_id = "X99")
  )
  common <- archive_common_series(a)
  expect_false("X99" %in% common$scenario_id)
  expect_true(all(c("A1", "A2") %in% common$scenario_id))
})

test_that("the trend plot builds without raising", {
  g <- plot_performance_trend(good_archive())
  expect_s3_class(g, "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(g)))
})

test_that("materiality is judged against scale_metric, not the metric itself", {
  # Bias sits near zero, so its own baseline is a useless denominator: a move
  # from -2e-04 to +1e-04 is a "157% change" and means nothing. Scaling by rmse
  # asks the question that matters -- is the shift large relative to the
  # estimator's own error?
  a <- good_archive(n_runs = 6, n_series = 2)
  a$rmse <- 0.05
  a$mean_bias <- stats::runif(nrow(a), -3e-4, 3e-4)
  latest <- max(a$run_date)
  # A large relative move in bias, but tiny against rmse.
  a$mean_bias[a$run_date == latest & a$scenario_id == "A1"] <- 6e-4

  self_scaled <- detect_regressions(a, metric = "mean_bias",
                                    scale_metric = "mean_bias")
  rmse_scaled <- detect_regressions(a, metric = "mean_bias",
                                    scale_metric = "rmse")
  expect_true(nrow(rmse_scaled) == 0)        # 0.0006 vs rmse 0.05 -- immaterial
  expect_true(nrow(self_scaled) >= nrow(rmse_scaled))

  # ...and a shift that IS large against rmse still surfaces.
  a$mean_bias[a$run_date == latest & a$scenario_id == "A1"] <- 0.03
  expect_true("A1" %in% detect_regressions(a, metric = "mean_bias",
                                           scale_metric = "rmse")$scenario_id)
})

test_that("a near-zero MAD cannot manufacture a regression", {
  # Cached targets are bit-identical between runs, so MAD is zero or float
  # noise. Dividing by it turned fourth-decimal changes into z-scores of
  # several hundred: 70 spurious rmse flags on the real archive.
  a <- good_archive(n_runs = 6, n_series = 2)
  a$rmse <- 0.02 + stats::runif(nrow(a), -1e-17, 1e-17)   # bit-identical-ish
  latest <- max(a$run_date)
  a$rmse[a$run_date == latest & a$scenario_id == "A1"] <- 0.0200001  # trivial
  expect_equal(nrow(detect_regressions(a)), 0L)
})
