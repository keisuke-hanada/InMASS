validation_row <- function(check, passed, failures = 0L, details = "") {
  data.frame(
    check = check,
    passed = as.logical(passed),
    failures = as.integer(failures),
    details = details,
    stringsAsFactors = FALSE
  )
}

finite_or_na <- function(x) {
  is.na(x) | is.finite(x)
}

make_diagnostic_summary <- function(results) {
  keys <- unique(results[c("scenario_id", "estimator", "formula_id")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- results$scenario_id == key$scenario_id &
      results$estimator == key$estimator &
      results$formula_id == key$formula_id
    d <- results[idx, , drop = FALSE]
    data.frame(
      key,
      n_replicates = nrow(d),
      n_nonconverged = sum(!d$converged),
      n_warning_rows = sum(!is.na(d$warnings) & nzchar(d$warnings)),
      warnings = paste(unique(d$warnings[!is.na(d$warnings) & nzchar(d$warnings)]), collapse = " || "),
      diagnostics = paste(unique(d$diagnostics[!is.na(d$diagnostics) & nzchar(d$diagnostics)]), collapse = " || "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_main_dgm_rules <- function(paths, nsim) {
  scenarios <- build_main_scenario_grid(nsim)
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    spec <- scenarios[i, , drop = FALSE]
    params_file <- file.path(scenario_raw_dir(paths, spec$scenario_id), "params.rds")
    if (!file.exists(params_file)) {
      return(data.frame(
        scenario_id = spec$scenario_id,
        external_sample_size_rule_ok = FALSE,
        covariate_dgm_ok = FALSE,
        outcome_dgm_ok = FALSE,
        details = "params.rds not found",
        stringsAsFactors = FALSE
      ))
    }
    params <- readRDS(params_file)
    n_external <- params$n_external
    external_ok <- length(n_external) == params$K &&
      all(n_external %% 2 == 0) &&
      all(n_external >= params$n) &&
      all(n_external <= 4 * params$n)
    covariate_ok <- identical(params$dist_x2, spec$covariate_distribution) &&
      params$dist_x2 %in% c("normal", "chi2")
    outcome_ok <- identical(as.numeric(params$beta), c(1, 2, -1, 0, 0.5)) &&
      identical(as.numeric(params$truth), 2)
    data.frame(
      scenario_id = spec$scenario_id,
      external_sample_size_rule_ok = external_ok,
      covariate_dgm_ok = covariate_ok,
      outcome_dgm_ok = outcome_ok,
      details = paste0(
        "External sizes are generated as 2*round(U(n/2, 2n)), ",
        "which gives even arm-balanced totals in [n, 4n]. ",
        "Normal x2 is N(mu_k, 1); chi-squared x2 is ((Z1^2+Z2^2)/2)+mu_k-1, ",
        "with mean mu_k and variance 1."
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_pilot_internal <- function(paths, nsim = 10L) {
  results_file <- file.path(paths$summary, sprintf("main_results_nsim%d.csv", nsim))
  summary_file <- file.path(paths$summary, sprintf("main_summary_nsim%d.csv", nsim))
  if (!file.exists(results_file) || !file.exists(summary_file)) {
    stop("Run the v2 pilot before internal validation.")
  }

  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  expected_scenarios <- build_main_scenario_grid(nsim)
  expected_estimators <- c("target_only", "inmass", "meta_regression", "plugin")
  expected_rows <- nrow(expected_scenarios) * nrow(main_analysis_formulas()) * length(expected_estimators) * nsim

  checks <- list()
  checks[[length(checks) + 1L]] <- validation_row(
    "expected_result_row_count",
    nrow(results) == expected_rows,
    abs(nrow(results) - expected_rows),
    sprintf("Observed %d rows; expected %d.", nrow(results), expected_rows)
  )

  converged <- results[results$converged, , drop = FALSE]
  finite_fail <- converged[
    !(is.finite(converged$estimate) &
        is.finite(converged$se) &
        is.finite(converged$ci_low) &
        is.finite(converged$ci_high) &
        converged$se >= 0),
    ,
    drop = FALSE
  ]
  checks[[length(checks) + 1L]] <- validation_row(
    "finite_estimable_quantities",
    nrow(finite_fail) == 0,
    nrow(finite_fail),
    "Converged rows must have finite estimate, SE, and CI values with non-negative SE."
  )

  ci_fail <- converged[
    !(converged$ci_low <= converged$estimate & converged$estimate <= converged$ci_high),
    ,
    drop = FALSE
  ]
  checks[[length(checks) + 1L]] <- validation_row(
    "ci_contains_estimate",
    nrow(ci_fail) == 0,
    nrow(ci_fail),
    "Requires ci_low <= estimate <= ci_high for converged rows."
  )

  nonconverged <- results[!results$converged, , drop = FALSE]
  na_fail <- nonconverged[
    !(is.na(nonconverged$estimate) &
        is.na(nonconverged$se) &
        is.na(nonconverged$ci_low) &
        is.na(nonconverged$ci_high)),
    ,
    drop = FALSE
  ]
  checks[[length(checks) + 1L]] <- validation_row(
    "nonconverged_quantities_are_na",
    nrow(na_fail) == 0,
    nrow(na_fail),
    "Rows marked non-converged should carry NA estimate, SE, and CI values."
  )

  summary_ok <- (
    finite_or_na(summary$bias) &
      finite_or_na(summary$mse) &
      finite_or_na(summary$power) &
      (is.na(summary$mse) | summary$mse >= 0) &
      (is.na(summary$power) | (summary$power >= 0 & summary$power <= 1)) &
      summary$n_converged <= summary$n_replicates
  )
  checks[[length(checks) + 1L]] <- validation_row(
    "summary_metrics_valid",
    all(summary_ok),
    sum(!summary_ok),
    "Bias must be finite/NA; MSE non-negative/NA; power in [0,1]/NA; n_converged <= n_replicates."
  )

  plugin_rows <- results[results$estimator == "plugin", , drop = FALSE]
  plugin_summary <- summary[summary$estimator == "plugin", , drop = FALSE]
  plugin_estimable <- plugin_rows[plugin_rows$allocation %in% c("1to1", "3to1"), , drop = FALSE]
  plugin_nonestimable <- plugin_rows[plugin_rows$allocation == "4to0", , drop = FALSE]
  plugin_diag_ok <- nrow(plugin_estimable) > 0 &&
    all(plugin_estimable$converged) &&
    all(grepl("mapping=arm-level predicted treated-control contrast", plugin_estimable$diagnostics, fixed = TRUE)) &&
    all(grepl("contrast=", plugin_estimable$diagnostics, fixed = TRUE)) &&
    all(grepl("beta0_hat=", plugin_estimable$diagnostics, fixed = TRUE)) &&
    all(grepl("beta_M_hat=", plugin_estimable$diagnostics, fixed = TRUE))
  checks[[length(checks) + 1L]] <- validation_row(
    "plugin_delta_method_contrast_recorded",
    plugin_diag_ok,
    if (plugin_diag_ok) 0L else sum(!grepl("contrast=", plugin_estimable$diagnostics, fixed = TRUE)),
    "Estimable plug-in rows must converge and record beta0_hat, beta_M_hat, and coefficient-name contrast diagnostics."
  )
  plugin_4to0_ok <- nrow(plugin_nonestimable) > 0 &&
    all(!plugin_nonestimable$converged) &&
    all(is.na(plugin_nonestimable$estimate)) &&
    all(is.na(plugin_nonestimable$se)) &&
    all(is.na(plugin_nonestimable$ci_low)) &&
    all(is.na(plugin_nonestimable$ci_high)) &&
    all(grepl("target control covariate summaries are unavailable", plugin_nonestimable$diagnostics, fixed = TRUE))
  checks[[length(checks) + 1L]] <- validation_row(
    "plugin_not_estimable_in_4to0",
    plugin_4to0_ok,
    if (plugin_4to0_ok) 0L else nrow(plugin_nonestimable),
    "Plug-in estimator must be explicit NA/non-converged under 4to0 because target control covariate summaries are unavailable."
  )
  plugin_1to1_3to1_ok <- nrow(plugin_estimable) > 0 &&
    all(plugin_estimable$converged) &&
    all(is.finite(plugin_estimable$estimate)) &&
    all(is.finite(plugin_estimable$se)) &&
    all(is.finite(plugin_estimable$ci_low)) &&
    all(is.finite(plugin_estimable$ci_high))
  checks[[length(checks) + 1L]] <- validation_row(
    "plugin_estimable_in_1to1_and_3to1",
    plugin_1to1_3to1_ok,
    if (plugin_1to1_3to1_ok) 0L else sum(!plugin_estimable$converged),
    "Plug-in estimator should remain estimable under 1to1 and 3to1 target designs."
  )
  plugin_power_ok <- nrow(plugin_summary) > 0 &&
    all(is.na(plugin_summary$power) | (plugin_summary$power >= 0 & plugin_summary$power <= 1))
  checks[[length(checks) + 1L]] <- validation_row(
    "plugin_power_from_ci_valid",
    plugin_power_ok,
    if (plugin_power_ok) 0L else sum(!(is.na(plugin_summary$power) | (plugin_summary$power >= 0 & plugin_summary$power <= 1))),
    "Plug-in empirical power is computed from the stored delta-method confidence interval and must lie in [0, 1]."
  )

  treatment_only <- results[
    results$allocation == "4to0" & results$estimator == "target_only",
    ,
    drop = FALSE
  ]
  target_only_ok <- nrow(treatment_only) > 0 &&
    all(!treatment_only$converged) &&
    all(is.na(treatment_only$estimate)) &&
    all(is.na(treatment_only$se)) &&
    all(is.na(treatment_only$ci_low)) &&
    all(is.na(treatment_only$ci_high))
  checks[[length(checks) + 1L]] <- validation_row(
    "target_only_not_estimable_in_4to0",
    target_only_ok,
    if (target_only_ok) 0L else nrow(treatment_only),
    "Treatment-only target design should have explicit NA target-only treatment estimates."
  )

  true_delta_ok <- "true_delta" %in% names(results) &&
    "true_delta" %in% names(summary) &&
    all(results$true_delta == 2) &&
    all(summary$true_delta == 2)
  checks[[length(checks) + 1L]] <- validation_row(
    "true_delta_is_two",
    true_delta_ok,
    if (true_delta_ok) 0L else 1L,
    "All v2 result and summary rows must carry true_delta = 2."
  )

  dgm_checks <- validate_main_dgm_rules(paths, nsim)
  dgm_ok <- all(dgm_checks$external_sample_size_rule_ok) &&
    all(dgm_checks$covariate_dgm_ok) &&
    all(dgm_checks$outcome_dgm_ok)
  checks[[length(checks) + 1L]] <- validation_row(
    "main_dgm_rules_documented_and_checked",
    dgm_ok,
    sum(!(dgm_checks$external_sample_size_rule_ok & dgm_checks$covariate_dgm_ok & dgm_checks$outcome_dgm_ok)),
    "Checks external sample-size rule, covariate DGM labels, outcome beta, and true delta."
  )

  checks <- do.call(rbind, checks)
  diagnostics <- make_diagnostic_summary(results)
  write_csv(checks, file.path(paths$summary, sprintf("pilot_internal_validation_nsim%d.csv", nsim)))
  write_csv(diagnostics, file.path(paths$summary, sprintf("pilot_diagnostics_nsim%d.csv", nsim)))
  write_csv(dgm_checks, file.path(paths$summary, sprintf("pilot_dgm_checks_nsim%d.csv", nsim)))

  list(checks = checks, diagnostics = diagnostics, dgm_checks = dgm_checks)
}
