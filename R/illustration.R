# One-dataset illustration for Figure 1, column 1.
#
# Simulates a single meta-analytic dataset of 50 DiD studies with a
# substantial common time trend, converts 25 of them to pre-post, and fits
# five models:
#
#   did_half_full   full model, the first half of the DiD studies alone
#   pp_half_naive   the second half as PP studies, analysed at face value
#                (time_trend = "fixed_zero"; requires allow_no_did = TRUE)
#   mixed_full   full model, 25 DiD + 25 PP
#   mixed_naive  naive model, 25 DiD + 25 PP
#   did_all_full   full model, all n_studies as DiD (the oracle)
#
# Returns tidy posterior draws of treatment_effect_mean for KDE panels,
# plus the true normalised effect. The seed is FIXED (set before inspecting
# any output) so this is one representative realisation, not a curated one.

library(metadid)

run_illustration <- function(seed = 495L, n_studies = 20L,
                             sigma_effect = 0.08, true_trend = -0.15,
                             pkg = NULL) {
  # `pkg` is a dependency-tracking sentinel, not data. Check the type first so
  # a mis-bound positional argument fails with a message that names the cause.
  stopifnot(
    "n_studies must be a single even number" =
      is.numeric(n_studies) && length(n_studies) == 1L && n_studies %% 2 == 0,
    "sigma_effect must be a single number" =
      is.numeric(sigma_effect) && length(sigma_effect) == 1L,
    "true_trend must be a single number" =
      is.numeric(true_trend) && length(true_trend) == 1L
  )
  half <- n_studies %/% 2L
  dgp <- list(
    n_studies     = n_studies,
    n_control     = 100L,
    n_treatment   = 100L,
    true_effect   = -0.15,
    sigma_effect  = sigma_effect,
    true_trend    = true_trend,
    sigma_trend   = 0.01,
    baseline_mean = 0.45,
    baseline_sd   = 0.02,
    rho           = 0.5
  )
  sim <- do.call(metadid::simulate_meta_did, c(dgp, list(seed = seed)))
  ids <- unique(sim$study_id)
  tp  <- attr(sim, "true_params")

  subset_sim <- function(keep) {
    s <- dplyr::filter(sim, study_id %in% keep)
    attr(s, "true_params") <- dplyr::filter(tp, study_id %in% keep)
    s
  }
  did_half <- metadid::as_summary_did(subset_sim(ids[seq_len(half)]))
  pp_half  <- metadid::as_summary_pp(subset_sim(ids[(half + 1L):n_studies]))
  all_did  <- metadid::as_summary_did(sim)
  mixed    <- dplyr::bind_rows(did_half, pp_half)

  opts <- c(MCMC_OPTS, list(seed = seed))
  fits <- list(
    did_half_full  = do.call(metadid::meta_did,
                          c(list(summary_data = did_half), opts)),
    pp_half_naive  = do.call(metadid::meta_did_general,
                          c(list(summary_data       = pp_half,
                                 time_trend         = "fixed_zero",
                                 baseline_imbalance = "fixed_zero",
                                 allow_no_did       = TRUE), opts)),
    mixed_full  = do.call(metadid::meta_did,
                          c(list(summary_data = mixed), opts)),
    mixed_naive = do.call(metadid::meta_did_general,
                          c(list(summary_data       = mixed,
                                 time_trend         = "fixed_zero",
                                 baseline_imbalance = "fixed_zero"), opts)),
    did_all_full  = do.call(metadid::meta_did,
                          c(list(summary_data = all_did), opts))
  )

  # True normalised effect on the estimand scale E[theta_i / b_i]; matches
  # build_true_params(): (mu_theta / mu_b) * (1 + CV_b^2), the second-order
  # Jensen term for between-study baseline variation (SI Appendix, S6).
  cv2 <- (dgp$baseline_sd / dgp$baseline_mean)^2
  truth <- (dgp$true_effect / dgp$baseline_mean) * (1 + cv2)

  purrr::imap_dfr(fits, function(f, nm) {
    tibble::tibble(
      model = nm,
      draw  = as.numeric(
        f$fit$draws("treatment_effect_mean", format = "draws_matrix"))
    )
  }) |>
    dplyr::mutate(truth = truth, seed = seed,
                  n_studies = n_studies, half = half,
                  sigma_effect = sigma_effect, true_trend = true_trend)
}
