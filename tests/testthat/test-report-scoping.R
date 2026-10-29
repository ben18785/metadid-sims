# The Executive Summary decides what counts as a defect. Two ways it can be
# wrong: counting scenarios the reader cannot inspect, and applying a
# calibration test to scenarios that were built to fail.

test_that("validation and figure-sweep categories partition the scenarios", {
  cats <- sort(unique(substr(names(SCENARIO_CONFIGS), 1L, 1L)))
  expect_equal(sort(c(VALIDATION_CATEGORIES, FIGURE_SWEEP_CATEGORIES)), cats)
  expect_equal(intersect(VALIDATION_CATEGORIES, FIGURE_SWEEP_CATEGORIES), character(0))
})

test_that("every category with a report section is a validation category", {
  # A category counted in the defect total must be one the reader can go and
  # look at. Figure sweeps are excluded from the count precisely because the
  # report has no section for them.
  qmd <- readLines(file.path(.sims_root, "reports/validation-report.qmd"), warn = FALSE)
  sections <- sub("^## ([A-Z])\\..*$", "\\1", grep("^## [A-Z]\\.", qmd, value = TRUE))
  expect_equal(setdiff(sections, VALIDATION_CATEGORIES), character(0))
})

test_that("contrast claims name scenarios and arms that actually exist", {
  reg <- scenario_contrasts()
  tbl <- build_scenario_settings_table()
  for (i in seq_len(nrow(reg))) {
    sid <- reg$scenario_id[i]
    expect_true(sid %in% names(SCENARIO_CONFIGS), info = sid)
    arms <- tbl$model_label[tbl$scenario_id == sid]
    for (a in c(reg$better[i], reg$worse[i])) {
      expect_true(a %in% arms,
                  info = paste0(sid, ": contrast names arm '", a,
                                "' but that scenario runs: ",
                                paste(arms, collapse = ", ")))
    }
    expect_false(identical(reg$better[i], reg$worse[i]), info = sid)
  }
})

test_that("check_contrasts evaluates the ordering it claims to", {
  agg <- tibble::tibble(
    scenario_id = c("X1", "X1"),
    parameter   = "treatment_effect_mean",
    model_label = c("zero_mean", "shared"),
    rmse        = c(0.02, 0.05)
  )
  out <- check_contrasts(agg)
  x1 <- out[out$scenario_id == "X1", ]
  expect_true(x1$holds)
  expect_equal(x1$better_rmse, 0.02)
  expect_equal(x1$worse_rmse, 0.05)

  # ...and reports a violation rather than quietly passing.
  agg$rmse <- c(0.09, 0.05)
  expect_false(check_contrasts(agg)$holds[check_contrasts(agg)$scenario_id == "X1"])
})

test_that("a missing arm yields NA, not a silent pass", {
  agg <- tibble::tibble(scenario_id = "X1", parameter = "treatment_effect_mean",
                        model_label = "zero_mean", rmse = 0.02)
  out <- check_contrasts(agg)
  expect_true(is.na(out$holds[out$scenario_id == "X1"]))
})

test_that("any open issue is stated, not an empty placeholder", {
  iss <- known_open_issues()
  # Zero rows is legitimate; a row with blank fields is not.
  expect_true(all(nzchar(iss$id) & nzchar(iss$affects) & nzchar(iss$summary)))
  expect_setequal(names(iss), c("id", "affects", "summary"))
})

test_that("the post-only residual bias is recorded as resolved, not hidden", {
  # It was traced to this repo's own generator: the post-only branch of
  # simulate_imbalance_grid omitted the time trend, so metadid's alpha + beta
  # -> alpha conversion was applied to data already on the alpha scale. The
  # guard that matters is unchanged -- a finding must not vanish into the
  # expectations registry, which is where real defects go to be forgotten.
  res <- resolved_issues()
  expect_true("post-only-residual-bias" %in% res$id)
  expect_true(all(nzchar(res$cause) & nzchar(res$fix)))
  expect_false("post-only-residual-bias" %in% known_open_issues()$id)
  reg <- scenario_expectations()
  expect_false(any(grepl("post-only residual", reg$reason, ignore.case = TRUE)))
})
