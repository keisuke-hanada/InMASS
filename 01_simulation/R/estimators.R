result_row <- function(scenario_id, replicate, estimator, formula_id, estimate = NA_real_,
                       se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                       converged = FALSE, diagnostics = "", warnings = "") {
  data.frame(
    scenario_id = scenario_id,
    replicate = as.integer(replicate),
    estimator = estimator,
    formula_id = formula_id,
    estimate = as.numeric(estimate),
    se = as.numeric(se),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    converged = as.logical(converged),
    diagnostics = diagnostics,
    warnings = warnings,
    stringsAsFactors = FALSE
  )
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = paste(unique(warnings), collapse = " | "))
}

collapse_messages <- function(...) {
  messages <- unlist(list(...), use.names = FALSE)
  messages <- messages[!is.na(messages) & nzchar(messages)]
  paste(unique(messages), collapse = " | ")
}

extract_regression_term <- function(fit, term = "x1k") {
  coef <- stats::coef(fit)
  idx <- match(term, names(coef))
  if (is.na(idx) || is.na(coef[idx])) {
    return(list(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  if (inherits(fit, "iwlm")) {
    ci <- stats::confint(fit)
    se <- sqrt(diag(fit$cov))[idx]
  } else {
    ci <- stats::confint(fit)
    coef_table <- summary(fit)$coefficients
    row_idx <- match(term, row.names(coef_table))
    if (is.na(row_idx)) {
      return(list(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
    }
    se <- coef_table[row_idx, 2]
  }
  list(estimate = coef[idx], se = se, ci_low = ci[idx, 1], ci_high = ci[idx, 2])
}

fit_ipd_replicate <- function(target_ipd, formula, scenario_id, replicate, formula_id, treatment_var = "x1k") {
  captured <- capture_warnings(
    tryCatch(stats::lm(stats::as.formula(formula), data = target_ipd), error = function(e) e)
  )
  fit <- captured$value
  if (inherits(fit, "error")) {
    return(result_row(
      scenario_id, replicate, "target_only", formula_id,
      diagnostics = conditionMessage(fit),
      warnings = captured$warnings
    ))
  }
  extracted <- capture_warnings(extract_regression_term(fit, term = treatment_var))
  ext <- extracted$value
  result_row(
    scenario_id, replicate, "target_only", formula_id,
    ext$estimate, ext$se, ext$ci_low, ext$ci_high,
    converged = !is.na(ext$estimate),
    diagnostics = ifelse(is.na(ext$estimate), paste(treatment_var, "coefficient unavailable"), ""),
    warnings = collapse_messages(captured$warnings, extracted$warnings)
  )
}

fit_meta_replicate <- function(target_ipd, data_mean, data_var, formula, scenario_id, replicate, formula_id, treatment_var = "x1k") {
  target_ad <- target_aggregate_rows(target_ipd, names(data_mean), replicate, arm_col = treatment_var)
  all_mean <- rbind(data_mean, target_ad[target_ad$var == "mean", , drop = FALSE])
  all_var <- rbind(data_var, target_ad[target_ad$var == "var", , drop = FALSE])
  captured <- capture_warnings(fit_meta_regression(all_mean, all_var, formula))
  fit <- captured$value
  ext <- extract_meta_coefficient(fit, coefficient = treatment_var)
  result_row(
    scenario_id, replicate, "meta_regression", formula_id,
    ext$estimate, ext$se, ext$ci_low, ext$ci_high,
    converged = isTRUE(fit$converged) && !is.na(ext$estimate),
    diagnostics = ifelse(isTRUE(fit$converged), "", fit$error %||% "meta-regression failed"),
    warnings = captured$warnings
  )
}

fit_plugin_replicate <- function(target_ipd, data_mean, data_var, formula, scenario_id, replicate, formula_id, allocation) {
  if (allocation == "4to0") {
    return(result_row(
      scenario_id, replicate, "plugin", formula_id,
      diagnostics = "plugin estimator not estimable for 4to0: target control covariate summaries are unavailable"
    ))
  }
  captured <- capture_warnings(fit_meta_regression(data_mean, data_var, formula))
  fit <- captured$value
  contrasted <- capture_warnings(plugin_contrast(fit, formula, target_ipd))
  plug <- contrasted$value
  result_row(
    scenario_id, replicate, "plugin", formula_id,
    plug$estimate, plug$se, plug$ci_low, plug$ci_high,
    converged = isTRUE(plug$converged),
    diagnostics = plug$diagnostics,
    warnings = collapse_messages(captured$warnings, contrasted$warnings)
  )
}

fit_inmass_replicate <- function(target_ipd, data_mean, data_var, formula, formula_ma,
                                 scenario_id, replicate, formula_id, allocation, strata,
                                 base_seed = 1234L, density_covariates = NULL,
                                 density_include_quadratic = TRUE) {
  ps_meta <- if (allocation == "1to1") 1 else 2
  seed <- replicate_seed(scenario_id, replicate, paste("inmass", formula_id), base_seed)
  captured <- capture_warnings(
    fit_inmass_core(
      formula, formula_ma, data_mean, data_var, target_ipd, strata, ps_meta, seed,
      density_covariates = density_covariates,
      density_include_quadratic = density_include_quadratic
    )
  )
  fit <- captured$value
  if (!isTRUE(fit$converged)) {
    return(result_row(
      scenario_id, replicate, "inmass", formula_id,
      diagnostics = fit$error,
      warnings = captured$warnings
    ))
  }
  treatment_var <- all.vars(stats::as.formula(formula))[2]
  extracted <- capture_warnings(extract_regression_term(fit$fit, term = treatment_var))
  ext <- extracted$value
  result_row(
    scenario_id, replicate, "inmass", formula_id,
    ext$estimate, ext$se, ext$ci_low, ext$ci_high,
    converged = !is.na(ext$estimate),
    diagnostics = paste(names(fit$diagnostics), fit$diagnostics, sep = "=", collapse = ";"),
    warnings = collapse_messages(captured$warnings, extracted$warnings)
  )
}

run_estimators_for_replicate <- function(replicate, dat, spec, formula_spec, base_seed = 1234L,
                                         treatment_var = NULL, density_covariates = NULL,
                                         density_include_quadratic = TRUE) {
  if (is.null(treatment_var)) treatment_var <- spec$treatment_var %||% "x1k"
  target_ipd <- dat$target_ipd[dat$target_ipd$nsim == replicate, , drop = FALSE]
  ad <- dat$strata_ad[dat$strata_ad$nsim == replicate, , drop = FALSE]
  data_mean <- ad[ad$var == "mean", , drop = FALSE]
  data_var <- ad[ad$var == "var", , drop = FALSE]

  rows <- list()
  rows[[1L]] <- fit_ipd_replicate(target_ipd, formula_spec$formula, spec$scenario_id, replicate, formula_spec$formula_id, treatment_var)
  rows[[2L]] <- fit_inmass_replicate(
    target_ipd, data_mean, data_var, formula_spec$formula, spec$formula_ma,
    spec$scenario_id, replicate, formula_spec$formula_id,
    spec$allocation, spec$K, base_seed,
    density_covariates = density_covariates,
    density_include_quadratic = density_include_quadratic
  )
  rows[[3L]] <- fit_meta_replicate(target_ipd, data_mean, data_var, formula_spec$formula, spec$scenario_id, replicate, formula_spec$formula_id, treatment_var)
  rows[[4L]] <- fit_plugin_replicate(
    target_ipd, data_mean, data_var, formula_spec$formula,
    spec$scenario_id, replicate, formula_spec$formula_id, spec$allocation
  )
  do.call(rbind, rows)
}

run_estimators_for_replicate_safe <- function(replicate, dat, spec, formula_spec, base_seed = 1234L,
                                             treatment_var = NULL, density_covariates = NULL,
                                             density_include_quadratic = TRUE) {
  tryCatch(
    run_estimators_for_replicate(
      replicate, dat, spec, formula_spec, base_seed,
      treatment_var = treatment_var,
      density_covariates = density_covariates,
      density_include_quadratic = density_include_quadratic
    ),
    error = function(e) {
      diagnostics <- paste("replicate-level estimator failure:", conditionMessage(e))
      do.call(rbind, lapply(
        c("target_only", "inmass", "meta_regression", "plugin"),
        function(estimator) result_row(
          spec$scenario_id, replicate, estimator, formula_spec$formula_id,
          diagnostics = diagnostics
        )
      ))
    }
  )
}

run_estimators_for_scenario <- function(dat, spec, formula_spec, base_seed = 1234L,
                                        n_workers = 1L, cluster = NULL) {
  nsim <- dat$params$nsim
  treatment_var <- spec$treatment_var %||% "x1k"
  density_covariates <- spec$density_covariates %||% NULL
  if (!is.null(density_covariates) && length(density_covariates) == 1L) {
    density_covariates <- strsplit(density_covariates, "+", fixed = TRUE)[[1]]
  }
  density_include_quadratic <- if (!is.null(spec$density_include_quadratic)) {
    isTRUE(spec$density_include_quadratic)
  } else {
    TRUE
  }
  rows <- parallel_lapply(
    seq_len(nsim),
    run_estimators_for_replicate_safe,
    dat = dat,
    spec = spec,
    formula_spec = formula_spec,
    base_seed = base_seed,
    treatment_var = treatment_var,
    density_covariates = density_covariates,
    density_include_quadratic = density_include_quadratic,
    n_workers = n_workers,
    cluster = cluster
  )
  do.call(rbind, rows)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
