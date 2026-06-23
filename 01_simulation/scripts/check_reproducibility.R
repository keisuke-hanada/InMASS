compare_output_files <- function(paths_a, paths_b) {
  files_a <- list.files(
    c(paths_a$raw, paths_a$summary),
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  root_a <- normalizePath(dirname(paths_a$raw), winslash = "/")
  root_b <- normalizePath(dirname(paths_b$raw), winslash = "/")
  norm_a <- gsub("\\\\", "/", normalizePath(files_a, winslash = "/"))
  rel_a <- substr(norm_a, nchar(root_a) + 2L, nchar(norm_a))

  files_b <- file.path(root_b, rel_a)
  rows <- lapply(seq_along(files_a), function(i) {
    file_a <- files_a[i]
    file_b <- files_b[i]
    exists_b <- file.exists(file_b)
    identical_object <- FALSE
    detail <- ""
    if (exists_b) {
      ext <- tolower(tools::file_ext(file_a))
      identical_object <- if (ext == "rds") {
        identical(readRDS(file_a), readRDS(file_b))
      } else if (ext == "csv") {
        identical(
          utils::read.csv(file_a, stringsAsFactors = FALSE, check.names = FALSE),
          utils::read.csv(file_b, stringsAsFactors = FALSE, check.names = FALSE)
        )
      } else {
        identical(readBin(file_a, "raw", file.info(file_a)$size),
                  readBin(file_b, "raw", file.info(file_b)$size))
      }
      if (!identical_object) detail <- "objects differ"
    } else {
      detail <- "missing in second run"
    }
    data.frame(
      relative_path = rel_a[i],
      exists_in_second_run = exists_b,
      identical = identical_object,
      detail = detail,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  files_b_all <- list.files(
    c(paths_b$raw, paths_b$summary),
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  norm_b <- gsub("\\\\", "/", normalizePath(files_b_all, winslash = "/"))
  rel_b_all <- substr(norm_b, nchar(root_b) + 2L, nchar(norm_b))
  extra_b <- setdiff(rel_b_all, rel_a)
  if (length(extra_b)) {
    out <- rbind(
      out,
      data.frame(
        relative_path = extra_b,
        exists_in_second_run = TRUE,
        identical = FALSE,
        detail = "extra in second run",
        stringsAsFactors = FALSE
      )
    )
  }
  out
}

check_pilot_reproducibility <- function(repo_root = getwd(), nsim = 2L, base_seed = 1234L,
                                        scenario_ids = NULL) {
  if (is.null(scenario_ids)) {
    scenario_ids <- c("main_1to1_K05_n020_normal", "main_4to0_K05_n020_chi2")
  }
  paths_a <- simulation_paths(repo_root, output_root = file.path("results", "repro_run_a"))
  paths_b <- simulation_paths(repo_root, output_root = file.path("results", "repro_run_b"))

  invisible(run_main_scenarios(paths_a, nsim = nsim, base_seed = base_seed, scenario_ids = scenario_ids))
  invisible(run_main_scenarios(paths_b, nsim = nsim, base_seed = base_seed, scenario_ids = scenario_ids))
  comparison <- compare_output_files(paths_a, paths_b)

  report_paths <- simulation_paths(repo_root)
  ensure_simulation_dirs(report_paths)
  write_csv(comparison, file.path(report_paths$summary, sprintf("pilot_reproducibility_nsim%d.csv", nsim)))
  comparison
}
