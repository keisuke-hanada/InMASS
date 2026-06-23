truncation_family_file <- function(paths, simulation_family, nsim) {
  prefix <- switch(
    simulation_family,
    main = "main",
    robustness_multicov = "robustness_multicov",
    robustness_nonlinear = "robustness_nonlinear",
    stop("Unknown simulation_family: ", simulation_family)
  )
  file.path(paths$summary, sprintf("%s_ripd_truncation_diagnostics_nsim%d.csv", prefix, nsim))
}

scenario_metadata_for_truncation <- function(spec, formula_spec, replicate) {
  list(
    simulation_family = spec$simulation_family %||% "main",
    dgm = spec$dgm %||% NA_character_,
    scenario_id = spec$scenario_id,
    replicate = as.integer(replicate),
    allocation = spec$allocation,
    K = as.integer(spec$K),
    n = as.integer(spec$n),
    covariate_distribution = spec$covariate_distribution %||% NA_character_,
    formula_id = formula_spec$formula_id
  )
}

collect_ripd_truncation_for_scenario <- function(dat, spec, formulas, base_seed = 1234L) {
  rows <- list()
  pos <- 1L
  for (j in seq_len(nrow(formulas))) {
    formula_spec <- as.list(formulas[j, , drop = FALSE])
    for (replicate in seq_len(dat$params$nsim)) {
      ad <- dat$strata_ad[dat$strata_ad$nsim == replicate, , drop = FALSE]
      data_mean <- ad[ad$var == "mean", , drop = FALSE]
      data_var <- ad[ad$var == "var", , drop = FALSE]
      meta_fit <- fit_meta_regression(data_mean, data_var, spec$formula_ma)
      seed <- replicate_seed(spec$scenario_id, replicate, paste("inmass", formula_spec$formula_id), base_seed)
      pseudo <- reconstruct_pseudo_ipd(
        data_mean,
        data_var,
        spec$formula_ma,
        meta_fit,
        as.integer(spec$K),
        seed,
        metadata = scenario_metadata_for_truncation(spec, formula_spec, replicate)
      )
      diag <- attr(pseudo, "truncation_diagnostics")
      if (!is.null(diag) && nrow(diag)) {
        rows[[pos]] <- diag
        pos <- pos + 1L
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

summarize_ripd_truncation <- function(diagnostics) {
  if (!nrow(diagnostics)) return(data.frame())
  diagnostics$replicate <- as.integer(diagnostics$replicate)
  diagnostics$K <- as.integer(diagnostics$K)
  diagnostics$n <- as.integer(diagnostics$n)
  diagnostics$truncated <- as.logical(diagnostics$truncated)
  keys <- c("simulation_family", "dgm", "allocation", "K", "n", "covariate_distribution", "formula_id")
  groups <- split(diagnostics, interaction(diagnostics[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(g) {
    data.frame(
      g[1, keys, drop = FALSE],
      n_replicates = length(unique(g$replicate)),
      n_reconstruction_units = nrow(g),
      n_truncated = sum(g$truncated, na.rm = TRUE),
      truncation_rate = mean(g$truncated, na.rm = TRUE),
      mean_untruncated_residual_variance = mean(g$resid_var_untruncated, na.rm = TRUE),
      min_untruncated_residual_variance = min(g$resid_var_untruncated, na.rm = TRUE),
      q05_untruncated_residual_variance = as.numeric(stats::quantile(g$resid_var_untruncated, 0.05, na.rm = TRUE, names = FALSE)),
      median_untruncated_residual_variance = stats::median(g$resid_var_untruncated, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

validate_ripd_truncation <- function(diagnostics, summary, expected_families = NULL) {
  checks <- list()
  checks[[length(checks) + 1L]] <- validation_row(
    "ripd_residual_variance_used_nonnegative",
    nrow(diagnostics) > 0 && all(diagnostics$resid_var_used >= 0),
    sum(diagnostics$resid_var_used < 0),
    "Requires resid_var_used >= 0 for every reconstruction unit."
  )
  indicator_ok <- diagnostics$truncated == (diagnostics$resid_var_untruncated <= 0)
  checks[[length(checks) + 1L]] <- validation_row(
    "ripd_truncation_indicator_consistent",
    nrow(diagnostics) > 0 && all(indicator_ok),
    sum(!indicator_ok),
    "Requires truncated iff resid_var_untruncated <= 0."
  )
  rate_ok <- summary$truncation_rate >= 0 & summary$truncation_rate <= 1
  checks[[length(checks) + 1L]] <- validation_row(
    "ripd_truncation_rate_in_unit_interval",
    nrow(summary) > 0 && all(rate_ok),
    sum(!rate_ok),
    "Summary truncation rates must be in [0, 1]."
  )
  positive_units <- summary$n_reconstruction_units > 0
  checks[[length(checks) + 1L]] <- validation_row(
    "ripd_positive_reconstruction_units",
    nrow(summary) > 0 && all(positive_units),
    sum(!positive_units),
    "Every summarized InMASS scenario/formula must have reconstruction units."
  )
  if (!is.null(expected_families)) {
    family_ok <- all(expected_families %in% unique(diagnostics$simulation_family))
    checks[[length(checks) + 1L]] <- validation_row(
      "ripd_expected_family_files_present",
      family_ok,
      sum(!expected_families %in% unique(diagnostics$simulation_family)),
      paste("Expected families:", paste(expected_families, collapse = ", "))
    )
  }
  do.call(rbind, checks)
}

write_ripd_truncation_outputs <- function(paths, nsim, family_diagnostics) {
  family_diagnostics <- family_diagnostics[vapply(family_diagnostics, nrow, integer(1)) > 0]
  diagnostics <- do.call(rbind, family_diagnostics)
  summary <- summarize_ripd_truncation(diagnostics)
  write_csv(summary, file.path(paths$summary, sprintf("ripd_truncation_summary_all_nsim%d.csv", nsim)))
  validation <- validate_ripd_truncation(
    diagnostics,
    summary,
    expected_families = names(family_diagnostics)
  )
  write_csv(validation, file.path(paths$summary, sprintf("ripd_truncation_validation_all_nsim%d.csv", nsim)))
  list(diagnostics = diagnostics, summary = summary, validation = validation)
}

summarize_existing_ripd_truncation <- function(paths, nsim = 10L) {
  families <- c("main", "robustness_multicov", "robustness_nonlinear")
  family_diagnostics <- list()
  for (family in families) {
    file <- truncation_family_file(paths, family, nsim)
    if (file.exists(file)) {
      family_diagnostics[[family]] <- utils::read.csv(file, stringsAsFactors = FALSE)
    }
  }
  if (!length(family_diagnostics)) {
    stop("No family-level RIPD truncation diagnostic files found under ", paths$summary)
  }
  write_ripd_truncation_outputs(paths, nsim, family_diagnostics)
}
