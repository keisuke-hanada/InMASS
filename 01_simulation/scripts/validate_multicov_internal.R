validate_multicov_aggregate_schema <- function(paths, nsim) {
  scenarios <- build_multicov_scenario_grid(nsim)
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    spec <- scenarios[i, , drop = FALSE]
    aggregate_file <- file.path(multicov_raw_dir(paths, spec$scenario_id), "aggregate_sample_replicate1.rds")
    params_file <- file.path(multicov_raw_dir(paths, spec$scenario_id), "params.rds")
    if (!file.exists(aggregate_file) || !file.exists(params_file)) {
      return(data.frame(
        scenario_id = spec$scenario_id,
        aggregate_schema_ok = FALSE,
        true_delta_ok = FALSE,
        details = "aggregate sample or params file missing",
        stringsAsFactors = FALSE
      ))
    }
    ad <- readRDS(aggregate_file)
    params <- readRDS(params_file)
    mean_rows <- ad[ad$var == "mean", , drop = FALSE]
    var_rows <- ad[ad$var == "var", , drop = FALSE]
    needed <- c("yik", "x1", "x2", "z", "strata", "nsim", "var", "n")
    schema_ok <- all(needed %in% names(ad)) &&
      nrow(mean_rows) > 0 &&
      nrow(var_rows) > 0 &&
      all(c("x1", "x2") %in% names(mean_rows)) &&
      all(c("x1", "x2") %in% names(var_rows))
    true_delta_ok <- identical(as.numeric(params$truth), 2) &&
      identical(params$dgm, "multicov") &&
      identical(params$simulation_family, "robustness_multicov")
    data.frame(
      scenario_id = spec$scenario_id,
      aggregate_schema_ok = schema_ok,
      true_delta_ok = true_delta_ok,
      details = "Aggregate rows include arm-specific means and variances for yik, x1, and x2.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_multicov_internal <- function(paths, nsim = 10L) {
  results_file <- file.path(paths$summary, sprintf("robustness_multicov_results_nsim%d.csv", nsim))
  summary_file <- file.path(paths$summary, sprintf("robustness_multicov_summary_nsim%d.csv", nsim))
  if (!file.exists(results_file) || !file.exists(summary_file)) {
    stop("Run the multicov pilot before internal validation.")
  }

  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  expected_rows <- nrow(build_multicov_scenario_grid(nsim)) * nrow(multicov_analysis_formulas()) * 4L * nsim
  checks <- list()
  checks[[length(checks) + 1L]] <- validation_row(
    "expected_multicov_result_row_count",
    nrow(results) == expected_rows,
    abs(nrow(results) - expected_rows),
    sprintf("Observed %d rows; expected %d.", nrow(results), expected_rows)
  )

  family_ok <- all(results$simulation_family == "robustness_multicov") &&
    all(summary$simulation_family == "robustness_multicov") &&
    all(results$dgm == "multicov") &&
    all(summary$dgm == "multicov")
  checks[[length(checks) + 1L]] <- validation_row(
    "simulation_family_and_dgm_are_separate",
    family_ok,
    if (family_ok) 0L else 1L,
    "All multicov outputs must be labeled robustness_multicov / multicov."
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
    "multicov_finite_estimable_quantities",
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
    "multicov_ci_contains_estimate",
    nrow(ci_fail) == 0,
    nrow(ci_fail),
    "Requires ci_low <= estimate <= ci_high for converged rows."
  )

  summary_ok <- (
    (is.na(summary$bias) | is.finite(summary$bias)) &
      (is.na(summary$mse) | (is.finite(summary$mse) & summary$mse >= 0)) &
      (is.na(summary$power) | (summary$power >= 0 & summary$power <= 1)) &
      summary$n_converged <= summary$n_replicates
  )
  checks[[length(checks) + 1L]] <- validation_row(
    "multicov_summary_metrics_valid",
    all(summary_ok),
    sum(!summary_ok),
    "Bias finite/NA; MSE non-negative/NA; power in [0,1]/NA."
  )

  target_4to0 <- results[results$allocation == "4to0" & results$estimator == "target_only", , drop = FALSE]
  plugin_4to0 <- results[results$allocation == "4to0" & results$estimator == "plugin", , drop = FALSE]
  plugin_estimable <- results[results$allocation %in% c("1to1", "3to1") & results$estimator == "plugin", , drop = FALSE]
  target_4to0_ok <- nrow(target_4to0) > 0 && all(!target_4to0$converged) && all(is.na(target_4to0$estimate))
  plugin_4to0_ok <- nrow(plugin_4to0) > 0 && all(!plugin_4to0$converged) && all(is.na(plugin_4to0$estimate)) &&
    all(grepl("target control covariate summaries are unavailable", plugin_4to0$diagnostics, fixed = TRUE))
  plugin_estimable_ok <- nrow(plugin_estimable) > 0 && all(plugin_estimable$converged)

  checks[[length(checks) + 1L]] <- validation_row("multicov_target_only_not_estimable_in_4to0", target_4to0_ok, if (target_4to0_ok) 0L else nrow(target_4to0), "")
  checks[[length(checks) + 1L]] <- validation_row("multicov_plugin_not_estimable_in_4to0", plugin_4to0_ok, if (plugin_4to0_ok) 0L else nrow(plugin_4to0), "")
  checks[[length(checks) + 1L]] <- validation_row("multicov_plugin_estimable_in_1to1_and_3to1", plugin_estimable_ok, if (plugin_estimable_ok) 0L else sum(!plugin_estimable$converged), "")

  true_delta_ok <- all(results$true_delta == 2) && all(summary$true_delta == 2)
  checks[[length(checks) + 1L]] <- validation_row("multicov_true_delta_is_two", true_delta_ok, if (true_delta_ok) 0L else 1L, "")

  schema <- validate_multicov_aggregate_schema(paths, nsim)
  schema_ok <- all(schema$aggregate_schema_ok) && all(schema$true_delta_ok)
  checks[[length(checks) + 1L]] <- validation_row(
    "multicov_aggregate_schema_and_dgm_valid",
    schema_ok,
    sum(!(schema$aggregate_schema_ok & schema$true_delta_ok)),
    "Checks aggregate x1/x2 means and variances plus DGM metadata."
  )

  checks <- do.call(rbind, checks)
  diagnostics <- make_diagnostic_summary(results)
  write_csv(checks, file.path(paths$summary, sprintf("robustness_multicov_internal_validation_nsim%d.csv", nsim)))
  write_csv(diagnostics, file.path(paths$summary, sprintf("robustness_multicov_diagnostics_nsim%d.csv", nsim)))
  write_csv(schema, file.path(paths$summary, sprintf("robustness_multicov_schema_checks_nsim%d.csv", nsim)))
  list(checks = checks, diagnostics = diagnostics, schema = schema)
}
