# Structural checks on SCENARIO_CONFIGS. These catch configuration mistakes at
# source time rather than ninety seconds into a fit -- or, worse, after a
# multi-hour pipeline run has produced results attributed to the wrong settings.

test_that("every scenario declares a resolvable simulator", {
  for (sid in names(SCENARIO_CONFIGS)) {
    dgp <- SCENARIO_CONFIGS[[sid]]$dgp
    expect_true(dgp$type %in% c("metadid", "bespoke"), info = sid)
    if (identical(dgp$type, "bespoke")) {
      expect_true(!is.null(dgp$bespoke_fn), info = sid)
      expect_true(is.function(tryCatch(match.fun(dgp$bespoke_fn),
                                       error = function(e) NULL)),
                  info = paste(sid, "-> unknown bespoke_fn:", dgp$bespoke_fn))
    }
  }
})

test_that("every comparison arm is labelled and labels are unique per scenario", {
  # An unlabelled arm silently becomes "default" and collides with another,
  # so two configurations get aggregated together under one name.
  for (sid in names(SCENARIO_CONFIGS)) {
    cmp <- SCENARIO_CONFIGS[[sid]]$compare
    if (is.null(cmp)) next
    labels <- vapply(cmp, function(a) a$label %||% NA_character_, character(1))
    expect_false(any(is.na(labels)), info = paste(sid, "has an unlabelled arm"))
    expect_equal(anyDuplicated(labels), 0L,
                 info = paste(sid, "has duplicate arm labels:",
                              paste(labels[duplicated(labels)], collapse = ", ")))
  }
})

test_that("every fit-config key is a real argument of the function it targets", {
  # This is the check that catches metadid's API drifting out from under the
  # sims. Passing an argument metadid no longer accepts fails here in
  # milliseconds instead of inside a fit; passing one it does not YET accept
  # (a stale installed build) fails here too.
  skip_if_not_installed("metadid")

  # Keys consumed by the harness itself rather than forwarded to metadid.
  harness_only <- c("label", "fn", "normalise", "data_format", "provide_rho",
                    "time_trend", "pp_likelihood", "priors",
                    "robust_heterogeneity", "design_effects",
                    "correlated_effects", "hierarchical_rho",
                    "covariates", "multiplicative_covariate")

  for (sid in names(SCENARIO_CONFIGS)) {
    cfg  <- SCENARIO_CONFIGS[[sid]]
    arms <- cfg$compare %||% list(list())
    for (arm in arms) {
      fit <- modifyList(cfg$fit, arm)
      fn  <- switch(fit$fn %||% "meta_did",
                    meta_did         = metadid::meta_did,
                    meta_did_general = metadid::meta_did_general,
                    stop("unknown fn in ", sid))
      unknown <- setdiff(names(fit), c(names(formals(fn)), harness_only))
      expect_equal(
        unknown, character(0),
        info = paste0(
          sid, " / ", arm$label %||% "default",
          ": fit config sets keys that ", fit$fn %||% "meta_did",
          "() does not accept: ", paste(unknown, collapse = ", "),
          ".\nIf these are arguments metadid gained recently, the INSTALLED ",
          "metadid is older than the one these scenarios were written for -- ",
          "reinstall it (R CMD INSTALL ../metadid), and in CI check that the ",
          "metadid branch being tested against matches this one. Installed ",
          "metadid: ", as.character(utils::packageVersion("metadid")), " at ",
          find.package("metadid"), "."
        )
      )
    }
  }
})

test_that("meta_did_general-only options are not handed to meta_did", {
  # time_trend and pp_likelihood are forwarded only for meta_did_general;
  # setting them on a meta_did arm silently does nothing.
  for (sid in names(SCENARIO_CONFIGS)) {
    cfg <- SCENARIO_CONFIGS[[sid]]
    for (arm in (cfg$compare %||% list(list()))) {
      fit <- modifyList(cfg$fit, arm)
      if (identical(fit$fn %||% "meta_did", "meta_did")) {
        expect_true(identical(fit$time_trend %||% "pooled", "pooled"),
                    info = paste(sid, arm$label %||% "default",
                                 "sets time_trend on a meta_did arm"))
      }
    }
  }
})

test_that("scenario_ids() reaches every scenario exactly once", {
  prefixes <- unique(substr(names(SCENARIO_CONFIGS), 1L, 1L))
  reached  <- unlist(lapply(prefixes, scenario_ids), use.names = FALSE)
  expect_equal(sort(reached), sort(names(SCENARIO_CONFIGS)))
  expect_equal(anyDuplicated(reached), 0L)
})

test_that("the settings table has one row per executed scenario-model pair", {
  tbl <- build_scenario_settings_table()
  expect_true(all(c("scenario_id", "model_label", "category") %in% names(tbl)))
  expect_equal(anyDuplicated(tbl[, c("scenario_id", "model_label")]), 0L)
  expect_setequal(unique(tbl$scenario_id), names(SCENARIO_CONFIGS))
})

test_that("expectation-registry entries point at scenarios that exist", {
  reg <- scenario_expectations()
  expect_equal(setdiff(reg$scenario_id, names(SCENARIO_CONFIGS)), character(0))
  # ...and at model labels those scenarios actually run, so an entry cannot
  # silently excuse nothing.
  tbl <- build_scenario_settings_table()
  labelled <- reg[!is.na(reg$model_label), ]
  for (i in seq_len(nrow(labelled))) {
    sid <- labelled$scenario_id[i]; lab <- labelled$model_label[i]
    expect_true(
      lab %in% tbl$model_label[tbl$scenario_id == sid],
      info = paste0("expectation for ", sid, " / ", lab,
                    " matches no model actually run for that scenario")
    )
  }
})

test_that("category X uses the randomisation-aware simulator throughout", {
  for (sid in scenario_ids("X")) {
    expect_equal(SCENARIO_CONFIGS[[sid]]$dgp$bespoke_fn,
                 "simulate_imbalance_grid", info = sid)
  }
})

test_that("the API-drift guard fires on an argument metadid does not accept", {
  # Confidence check on the guard itself. This is the failure that actually
  # occurred: the sims passed allow_unidentified_kappa while the INSTALLED
  # metadid predated it, and nothing noticed until a fit died mid-run with
  # "unused argument". The guard must catch that at source time.
  skip_if_not_installed("metadid")
  harness_only <- c("label", "fn", "normalise", "data_format", "provide_rho",
                    "time_trend", "pp_likelihood", "priors",
                    "robust_heterogeneity", "design_effects",
                    "correlated_effects", "hierarchical_rho",
                    "covariates", "multiplicative_covariate")
  bogus   <- modifyList(default_fit, list(a_setting_metadid_never_had = TRUE))
  unknown <- setdiff(names(bogus),
                     c(names(formals(metadid::meta_did)), harness_only))
  expect_equal(unknown, "a_setting_metadid_never_had")

  # ...and passes cleanly for the real default config.
  ok <- setdiff(names(default_fit),
                c(names(formals(metadid::meta_did)), harness_only))
  expect_equal(ok, character(0))
})
