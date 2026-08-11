# One-dataset illustration for Figure 1, column 1.
#
# Simulates a single meta-analytic dataset of 50 DiD studies with a
# substantial common time trend, converts 25 of them to pre-post, and fits
# five models:
#
#   did25_full   full model, the 25 remaining DiD studies alone
#   pp25_naive   the 25 PP studies alone, analysed at face value
#                (time_trend = "fixed_zero"; requires allow_no_did = TRUE)
#   mixed_full   full model, 25 DiD + 25 PP
#   mixed_naive  naive model, 25 DiD + 25 PP
#   did50_full   full model, all 50 original DiD studies (the oracle)
#
# Returns tidy posterior draws of treatment_effect_mean for KDE panels,
# plus the true normalised effect. The seed is FIXED (set before inspecting
# any output) so this is one representative realisation, not a curated one.

library(metadid)

run_illustration <- function(seed = 495L, pkg = NULL) {
  dgp <- list(
    n_studies     = 50L,
    n_control     = 100L,
    n_treatment   = 100L,
    true_effect   = -0.15,
    sigma_effect  = 0.03,
    true_trend    = -0.10,
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
  did_half <- metadid::as_summary_did(subset_sim(ids[1:25]))
  pp_half  <- metadid::as_summary_pp(subset_sim(ids[26:50]))
  all_did  <- metadid::as_summary_did(sim)
  mixed    <- dplyr::bind_rows(did_half, pp_half)

  opts <- c(MCMC_OPTS, list(seed = seed))
  fits <- list(
    did25_full  = do.call(metadid::meta_did,
                          c(list(summary_data = did_half), opts)),
    pp25_naive  = do.call(metadid::meta_did_general,
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
    did50_full  = do.call(metadid::meta_did,
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
    dplyr::mutate(truth = truth, seed = seed)
}
