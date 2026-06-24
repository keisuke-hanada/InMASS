run_ripd_truncation_for_family <- function(paths, simulation_family, nsim = 10L, base_seed = 1234L,
                                           n_workers = 1L) {
  ensure_simulation_dirs(paths)
  started_at <- Sys.time()
  message(sprintf("Starting RIPD truncation diagnostics for %s: nsim=%d, n_workers=%d", simulation_family, nsim, n_workers))
  cluster <- make_parallel_cluster(n_workers)
  on.exit(stop_parallel_cluster(cluster), add = TRUE)

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
    dat <- generator(spec_list, scenario_seed(spec$scenario_id, "data", base_seed), n_workers = n_workers, cluster = cluster)
    rows[[i]] <- collect_ripd_truncation_for_scenario(dat, spec_list, formulas, base_seed, n_workers = n_workers, cluster = cluster)
  }

  diagnostics <- do.call(rbind, rows)
  write_csv(diagnostics, truncation_family_file(paths, simulation_family, nsim))
  finished_at <- Sys.time()
  append_parallel_timing(paths, paste0(simulation_family, "_ripd_truncation"), nsim, n_workers, started_at, finished_at)
  message(sprintf("Finished RIPD truncation diagnostics for %s in %.1f seconds", simulation_family, as.numeric(difftime(finished_at, started_at, units = "secs"))))
  diagnostics
}

run_ripd_truncation_diagnostics <- function(paths, nsim = 10L, base_seed = 1234L,
                                            families = c("main", "robustness_multicov", "robustness_nonlinear"),
                                            n_workers = 1L) {
  family_diagnostics <- list()
  for (family in families) {
    family_diagnostics[[family]] <- run_ripd_truncation_for_family(paths, family, nsim, base_seed, n_workers = n_workers)
  }
  write_ripd_truncation_outputs(paths, nsim, family_diagnostics)
}
