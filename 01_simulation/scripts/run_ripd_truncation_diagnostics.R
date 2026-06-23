run_ripd_truncation_for_family <- function(paths, simulation_family, nsim = 10L, base_seed = 1234L) {
  ensure_simulation_dirs(paths)

  if (simulation_family == "main") {
    scenarios <- build_main_scenario_grid(nsim)
    formulas <- main_analysis_formulas()
    generator <- generate_main_data
  } else if (simulation_family == "robustness_multicov") {
    scenarios <- build_multicov_scenario_grid(nsim)
    formulas <- multicov_analysis_formulas()
    generator <- generate_multicov_data
  } else if (simulation_family == "robustness_nonlinear") {
    scenarios <- build_nonlinear_scenario_grid(nsim)
    formulas <- nonlinear_analysis_formulas()
    generator <- generate_nonlinear_data
  } else {
    stop("Unknown simulation_family: ", simulation_family)
  }

  rows <- list()
  for (i in seq_len(nrow(scenarios))) {
    spec <- scenarios[i, , drop = FALSE]
    spec_list <- as.list(spec)
    message(sprintf("Diagnosing RIPD truncation for %s (%d/%d)", spec$scenario_id, i, nrow(scenarios)))
    dat <- generator(spec_list, scenario_seed(spec$scenario_id, "data", base_seed))
    rows[[i]] <- collect_ripd_truncation_for_scenario(dat, spec_list, formulas, base_seed)
  }

  diagnostics <- do.call(rbind, rows)
  write_csv(diagnostics, truncation_family_file(paths, simulation_family, nsim))
  diagnostics
}

run_ripd_truncation_diagnostics <- function(paths, nsim = 10L, base_seed = 1234L,
                                            families = c("main", "robustness_multicov", "robustness_nonlinear")) {
  family_diagnostics <- list()
  for (family in families) {
    family_diagnostics[[family]] <- run_ripd_truncation_for_family(paths, family, nsim, base_seed)
  }
  write_ripd_truncation_outputs(paths, nsim, family_diagnostics)
}
