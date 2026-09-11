# The failure mode this file exists for: a parameter the model reports that
# NOTHING ever scores.
#
# Before the randomisation work the entire baseline-imbalance block --
# baseline_difference_mean, baseline_difference_sd -- was absent from
# assess_one()'s param_map. The P-category scenarios measured gamma's knock-on
# effect on the pooled treatment effect but never checked whether gamma itself
# was recovered, and nothing anywhere said so. A silently unscored parameter
# looks exactly like a scored one that always passes.
#
# The guard: every population parameter the fitted model reports must either be
# scored, or appear on an explicit list with a stated reason.

# assess_one() returns a COLUMN-LESS tibble when nothing scores, so `$parameter`
# warns rather than returning character(0). Access it safely.
scored_names <- function(res) {
  if (!"parameter" %in% names(res)) character(0) else as.character(res$parameter)
}

# Population-level parameters metadid can report, per R/fit.R's extract_posteriors()
# and metadid's summary() method.
REPORTED <- c(
  "treatment_effect_mean", "treatment_effect_sd",
  "time_trend_mean", "time_trend_sd",
  "baseline_control_mean",
  "baseline_difference_mean", "baseline_difference_sd", "kappa"
)

# Deliberately unscored, with reasons. Adding to this list should require an
# argument, which is the point of making it explicit.
UNSCORED_BY_DESIGN <- c(
  # Reported on the raw scale only; scored via baseline_mean when
  # normalise = FALSE, and not estimated at all when normalise = TRUE.
  baseline_control_mean =
    "scored as baseline_mean in the unnormalised branch; not a free parameter when normalised"
)

test_that("every reported population parameter is scored or explicitly excused", {
  mk_post <- function(nm) tibble(parameter = nm, mean = 0.5, sd = 0.1,
                                 lo = 0.4, hi = 0.6, max_rhat = 1,
                                 min_ess_bulk = 1000, n_divergent = 0L,
                                 n_max_treedepth = 0L, min_ebfmi = 1)
  posteriors <- bind_rows(lapply(REPORTED, mk_post))

  # A truth table carrying every column build_true_params() can emit, so the
  # only reason a parameter goes unscored is that param_map omits it.
  truths <- tibble(
    treatment_effect_mean_normalised = 0.5, treatment_effect_sd_normalised = 0.5,
    time_trend_mean_normalised = 0.5,       time_trend_sd_normalised = 0.5,
    treatment_effect_mean_raw = 0.5,        treatment_effect_sd_raw = 0.5,
    time_trend_mean_raw = 0.5,              time_trend_sd_raw = 0.5,
    baseline_difference_mean_normalised = 0.5, baseline_difference_sd_normalised = 0.5,
    baseline_difference_mean_raw = 0.5,        baseline_difference_sd_raw = 0.5,
    baseline_mean = 0.45, kappa = 0.5
  )

  for (norm in c(TRUE, FALSE)) {
    scored <- scored_names(assess_one(posteriors, truths,
                           list(normalise = norm, mu_gamma = "estimated")))
    missed <- setdiff(REPORTED, c(scored, names(UNSCORED_BY_DESIGN)))
    expect_equal(
      missed, character(0),
      info = paste0(
        "normalise = ", norm, ": these parameters are reported by the model but ",
        "scored by nothing, and are not on the explicit exclusion list: ",
        paste(missed, collapse = ", "),
        ". Either add them to assess_one()'s param_map or add them to ",
        "UNSCORED_BY_DESIGN with a reason."
      )
    )
  }
})

test_that("the exclusion list has no stale entries", {
  # An excused parameter that the model no longer reports means the list has
  # drifted and is hiding nothing.
  expect_equal(setdiff(names(UNSCORED_BY_DESIGN), REPORTED), character(0))
})

test_that("the guard actually fires when a parameter goes unscored", {
  # Confidence check on the guard itself: pretend the model reports something
  # assess_one() knows nothing about, and confirm it is detected.
  post <- tibble(parameter = "a_parameter_nothing_scores", mean = 1, sd = 0.1,
                 lo = 0.9, hi = 1.1, max_rhat = 1, min_ess_bulk = 1000,
                 n_divergent = 0L, n_max_treedepth = 0L, min_ebfmi = 1)
  truths <- tibble(treatment_effect_mean_normalised = 0.5,
                   treatment_effect_sd_normalised = 0.5,
                   time_trend_mean_normalised = 0.5,
                   time_trend_sd_normalised = 0.5, baseline_mean = 0.45)
  scored <- scored_names(assess_one(post, truths, list(normalise = TRUE)))
  expect_false("a_parameter_nothing_scores" %in% scored)
  expect_equal(scored, character(0))
})
