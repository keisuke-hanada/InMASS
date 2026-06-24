run_main_scenarios <- function(paths, nsim = 10L, base_seed = 1234L, scenario_ids = NULL,
                               n_workers = 1L, include_ripd_truncation = FALSE) {
  ensure_simulation_dirs(paths)
  started_at <- Sys.time()
  message(sprintf("Starting main scenarios: nsim=%d, n_workers=%d", nsim, n_workers))
  cluster <- make_parallel_cluster(n_workers)
  on.exit(stop_parallel_cluster(cluster), add = TRUE)
  scenarios <- build_main_scenario_grid(nsim)
  formulas <- main_analysis_formulas()
  if (!is.null(scenario_ids)) {
    scenarios <- scenarios[scenarios$scenario_id %in% scenario_ids, , drop = FALSE]
  }

  all_results <- list()
  all_summaries <- list()
  all_truncation <- if (include_ripd_truncation) list() else NULL
  result_pos <- 1L
  summary_pos <- 1L

  for (i in seq_len(nrow(scenarios))) {
    spec <- scenarios[i, , drop = FALSE]
    spec_list <- as.list(spec)
    message(sprintf("Running %s (%d/%d)", spec$scenario_id, i, nrow(scenarios)))

    dat <- generate_main_data(spec_list, scenario_seed(spec$scenario_id, "data", base_seed), n_workers = n_workers, cluster = cluster)
    save_rds(dat$params, file.path(scenario_raw_dir(paths, spec$scenario_id), "params.rds"))

    scenario_results <- list()
    for (j in seq_len(nrow(formulas))) {
      formula_spec <- as.list(formulas[j, , drop = FALSE])
      res <- run_estimators_for_scenario(dat, spec_list, formula_spec, base_seed, n_workers = n_workers, cluster = cluster)
      scenario_results[[j]] <- res
      save_rds(
        res,
        file.path(scenario_raw_dir(paths, spec$scenario_id), paste0("results_", formula_spec$formula_id, ".rds"))
      )
    }
    if (include_ripd_truncation) {
      all_truncation[[i]] <- collect_ripd_truncation_for_scenario(dat, spec_list, formulas, base_seed, n_workers = n_workers, cluster = cluster)
    }

    scenario_results <- do.call(rbind, scenario_results)
    scenario_results$true_delta <- spec$truth
    scenario_results$allocation <- spec$allocation
    scenario_results$K <- spec$K
    scenario_results$n <- spec$n
    scenario_results$covariate_distribution <- spec$covariate_distribution
    all_results[[result_pos]] <- scenario_results
    result_pos <- result_pos + 1L

    scenario_summary <- evaluate_estimates(scenario_results, truth = spec$truth)
    scenario_summary$allocation <- spec$allocation
    scenario_summary$K <- spec$K
    scenario_summary$n <- spec$n
    scenario_summary$covariate_distribution <- spec$covariate_distribution
    scenario_summary$true_delta <- spec$truth
    all_summaries[[summary_pos]] <- scenario_summary
    summary_pos <- summary_pos + 1L
    write_csv(scenario_summary, file.path(scenario_raw_dir(paths, spec$scenario_id), "summary.csv"))
  }

  results <- do.call(rbind, all_results)
  summary <- do.call(rbind, all_summaries)
  write_csv(results, file.path(paths$summary, sprintf("main_results_nsim%d.csv", nsim)))
  write_csv(summary, file.path(paths$summary, sprintf("main_summary_nsim%d.csv", nsim)))
  if (include_ripd_truncation) {
    truncation <- do.call(rbind, all_truncation)
    write_csv(truncation, truncation_family_file(paths, "main", nsim))
  }
  finished_at <- Sys.time()
  append_parallel_timing(paths, "main", nsim, n_workers, started_at, finished_at)
  message(sprintf("Finished main scenarios in %.1f seconds", as.numeric(difftime(finished_at, started_at, units = "secs"))))
  list(results = results, summary = summary)
}
