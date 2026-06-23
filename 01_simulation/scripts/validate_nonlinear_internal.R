validate_nonlinear_aggregate_schema <- function(paths, nsim) {
  scenarios <- build_nonlinear_scenario_grid(nsim)
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    spec <- scenarios[i, , drop = FALSE]
    aggregate_file <- file.path(nonlinear_raw_dir(paths, spec$scenario_id), "aggregate_sample_replicate1.rds")
    params_file <- file.path(nonlinear_raw_dir(paths, spec$scenario_id), "params.rds")
    if (!file.exists(aggregate_file) || !file.exists(params_file)) {
      return(data.frame(
        scenario_id = spec$scenario_id,
        aggregate_schema_ok = FALSE,
        moment_derivation_ok = FALSE,
        true_delta_ok = FALSE,
        details = "aggregate sample or params file missing",
        stringsAsFactors = FALSE
      ))
    }

    ad <- readRDS(aggregate_file)
    params <- readRDS(params_file)
    mean_rows <- ad[ad$var == "mean", , drop = FALSE]
    var_rows <- ad[ad$var == "var", , drop = FALSE]
    needed <- c("yik", "x", "z", "strata", "nsim", "var", "n", "x_second", "x_centered_second")
    schema_ok <- all(needed %in% names(ad)) && nrow(mean_rows) > 0 && nrow(var_rows) > 0

    moment_ok <- FALSE
    if (schema_ok) {
      key_mean <- paste(mean_rows$strata, mean_rows$z, mean_rows$nsim, sep = "::")
      key_var <- paste(var_rows$strata, var_rows$z, var_rows$nsim, sep = "::")
      x_var <- var_rows$x
      names(x_var) <- key_var
      expected_second <- x_var[key_mean] + mean_rows$x^2
      expected_centered <- expected_second - 1
      moment_ok <- all(is.finite(mean_rows$x_second)) &&
        all(is.finite(mean_rows$x_centered_second)) &&
        all(abs(mean_rows$x_second - as.numeric(expected_second)) < 1e-10) &&
        all(abs(mean_rows$x_centered_second - as.numeric(expected_centered)) < 1e-10) &&
        all(is.na(var_rows$x_second)) &&
        all(is.na(var_rows$x_centered_second)) &&
        isTRUE(params$nonlinear_derived_features_from_aggregate)
    }

    true_delta_ok <- identical(as.numeric(params$truth), 2) &&
      identical(params$dgm, "nonlinear") &&
      identical(params$simulation_family, "robustness_nonlinear")
    data.frame(
      scenario_id = spec$scenario_id,
      aggregate_schema_ok = schema_ok,
      moment_derivation_ok = moment_ok,
      true_delta_ok = true_delta_ok,
      details = "Aggregate rows include x mean/variance and derive E[X^2] and E[X^2-1] from them.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_nonlinear_internal <- function(paths, nsim = 10L) {
  results_file <- file.path(paths$summary, sprintf("robustness_nonlinear_results_nsim%d.csv", nsim))
  summary_file <- file.path(paths$summary, sprintf("robustness_nonlinear_summary_nsim%d.csv", nsim))
  if (!file.exists(results_file) || !file.exists(summary_file)) {
    stop("Run the nonlinear pilot before internal validation.")
  }

  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  expected_rows <- nrow(build_nonlinear_scenario_grid(nsim)) * nrow(nonlinear_analysis_formulas()) * 4L * nsim
  checks <- list()
  checks[[length(checks) + 1L]] <- validation_row(
    "expected_nonlinear_result_row_count",
    nrow(results) == expected_rows,
    abs(nrow(results) - expected_rows),
    sprintf("Observed %d rows; expected %d.", nrow(results), expected_rows)
  )

  family_ok <- all(results$simulation_family == "robustness_nonlinear") &&
    all(summary$simulation_family == "robustness_nonlinear") &&
    all(results$dgm == "nonlinear") &&
    all(summary$dgm == "nonlinear")
  checks[[length(checks) + 1L]] <- validation_row(
    "simulation_family_and_dgm_are_separate",
    family_ok,
    if (family_ok) 0L else 1L,
    "All nonlinear outputs must be labeled robustness_nonlinear / nonlinear."
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
    "nonlinear_finite_estimable_quantities",
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
    "nonlinear_ci_contains_estimate",
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
    "nonlinear_summary_metrics_valid",
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

  checks[[length(checks) + 1L]] <- validation_row("nonlinear_target_only_not_estimable_in_4to0", target_4to0_ok, if (target_4to0_ok) 0L else nrow(target_4to0), "")
  checks[[length(checks) + 1L]] <- validation_row("nonlinear_plugin_not_estimable_in_4to0", plugin_4to0_ok, if (plugin_4to0_ok) 0L else nrow(plugin_4to0), "")
  checks[[length(checks) + 1L]] <- validation_row("nonlinear_plugin_estimable_in_1to1_and_3to1", plugin_estimable_ok, if (plugin_estimable_ok) 0L else sum(!plugin_estimable$converged), "")

  true_delta_ok <- all(results$true_delta == 2) && all(summary$true_delta == 2)
  checks[[length(checks) + 1L]] <- validation_row("nonlinear_true_delta_is_two", true_delta_ok, if (true_delta_ok) 0L else 1L, "")

  schema <- validate_nonlinear_aggregate_schema(paths, nsim)
  schema_ok <- all(schema$aggregate_schema_ok) && all(schema$moment_derivation_ok) && all(schema$true_delta_ok)
  checks[[length(checks) + 1L]] <- validation_row(
    "nonlinear_aggregate_schema_moments_and_dgm_valid",
    schema_ok,
    sum(!(schema$aggregate_schema_ok & schema$moment_derivation_ok & schema$true_delta_ok)),
    "Checks aggregate x mean/variance, derived nonlinear moments, and DGM metadata."
  )

  d <- prepare_legacy_plot_data(results)
  d$model <- "nonlinear"
  power_data <- aggregate_power_curve_legacy(d, "Method_main", b.max = 4, h = 0.1)
  power_plot_data <- power_data[power_data$strata == 10, , drop = FALSE]
  power_diag <- validate_power_plot_data(
    power_plot_data,
    key_cols = c("n", "formula", "Method", "sd"),
    path = file.path(paths$summary, sprintf("robustness_nonlinear_power_plot_diagnostics_nsim%d.csv", nsim))
  )
  checks[[length(checks) + 1L]] <- validation_row(
    "nonlinear_power_plot_data_valid",
    all(power_diag$passed),
    sum(!power_diag$passed),
    "Power plotting data must be summarized and contain no replicate column or duplicate keys."
  )

  checks <- do.call(rbind, checks)
  diagnostics <- make_diagnostic_summary(results)
  write_csv(checks, file.path(paths$summary, sprintf("robustness_nonlinear_internal_validation_nsim%d.csv", nsim)))
  write_csv(diagnostics, file.path(paths$summary, sprintf("robustness_nonlinear_diagnostics_nsim%d.csv", nsim)))
  write_csv(schema, file.path(paths$summary, sprintf("robustness_nonlinear_schema_checks_nsim%d.csv", nsim)))
  list(checks = checks, diagnostics = diagnostics, schema = schema, power_diagnostics = power_diag)
}
