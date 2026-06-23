simulation_paths <- function(repo_root = getwd(), output_root = "results") {
  sim_root <- file.path(repo_root, "01_simulation")
  results_root <- file.path(repo_root, output_root)
  list(
    repo_root = repo_root,
    sim_root = sim_root,
    raw = file.path(results_root, "raw"),
    summary = file.path(results_root, "summary"),
    figures_main = file.path(repo_root, "figures", "main"),
    figures_supplement = file.path(repo_root, "figures", "supplement"),
    v1_root = file.path(repo_root, "04_simulation-v1.0")
  )
}

ensure_simulation_dirs <- function(paths) {
  dirs <- c(paths$raw, paths$summary, paths$figures_main, paths$figures_supplement)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}
