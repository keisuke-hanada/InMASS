sort_estimator_output <- function(x) {
  key <- intersect(
    c("scenario_id", "replicate", "estimator", "formula_id", "allocation", "K", "n", "covariate_distribution"),
    names(x)
  )
  x[do.call(order, x[key]), , drop = FALSE]
}

compare_parallel_outputs <- function(seq_file, par_file, tolerance = 1e-10) {
  seq_out <- sort_estimator_output(utils::read.csv(seq_file, stringsAsFactors = FALSE))
  par_out <- sort_estimator_output(utils::read.csv(par_file, stringsAsFactors = FALSE))
  common <- intersect(names(seq_out), names(par_out))
  seq_out <- seq_out[common]
  par_out <- par_out[common]
  identical_dims <- identical(dim(seq_out), dim(par_out))
  numeric_cols <- names(seq_out)[vapply(seq_out, is.numeric, logical(1))]
  char_cols <- setdiff(names(seq_out), numeric_cols)
  numeric_ok <- all(vapply(numeric_cols, function(col) {
    isTRUE(all.equal(seq_out[[col]], par_out[[col]], tolerance = tolerance, check.attributes = FALSE))
  }, logical(1)))
  nonnumeric_ok <- identical(seq_out[char_cols], par_out[char_cols])
  data.frame(
    check = c("identical_dimensions", "numeric_columns_within_tolerance", "nonnumeric_columns_identical"),
    passed = c(identical_dims, numeric_ok, nonnumeric_ok),
    failures = c(
      if (identical_dims) 0L else 1L,
      if (numeric_ok) 0L else 1L,
      if (nonnumeric_ok) 0L else 1L
    ),
    stringsAsFactors = FALSE
  )
}

check_parallel_reproducibility <- function(repo_root = getwd(), nsim = 3L, base_seed = 1234L) {
  seq_paths <- simulation_paths(repo_root, output_root = "results_parallel_repro_seq")
  par_paths <- simulation_paths(repo_root, output_root = "results_parallel_repro_par")
  scenario_ids <- "main_1to1_K05_n020_normal"
  run_main_scenarios(seq_paths, nsim = nsim, base_seed = base_seed, scenario_ids = scenario_ids, n_workers = 1L)
  run_main_scenarios(par_paths, nsim = nsim, base_seed = base_seed, scenario_ids = scenario_ids, n_workers = 2L)
  checks <- compare_parallel_outputs(
    file.path(seq_paths$summary, sprintf("main_results_nsim%d.csv", nsim)),
    file.path(par_paths$summary, sprintf("main_results_nsim%d.csv", nsim))
  )
  write_csv(checks, file.path(par_paths$summary, sprintf("parallel_reproducibility_check_nsim%d.csv", nsim)))
  checks
}
