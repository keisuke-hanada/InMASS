# Final simulation workflow wrapper.
#
# Run from the repository root with:
#   Rscript 01_simulation/run_final_simulations.R
#
# The script also works when run from inside 01_simulation/ as:
#   Rscript run_final_simulations.R
#
# Defaults are intentionally small for testing. For the final manuscript run,
# change NSIM to 10000. Adjust N_WORKERS if memory pressure or worker failures
# occur on the target machine.

NSIM <- 10
N_WORKERS <- 20

OUTPUT_ROOT_MAIN <- "results_final_main"
OUTPUT_ROOT_MULTICOV <- "results_final_multicov"
OUTPUT_ROOT_NONLINEAR <- "results_final_nonlinear"
OUTPUT_ROOT_TRUNC <- "results_final_trunc"
LOG_DIR <- "logs_final"

# Output roots are used exactly as written above, including when NSIM = 10.
# This keeps test and final runs reproducible: change the constants if a
# different output naming convention is desired.

detect_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path("01_simulation", "run_final_simulations.R"), winslash = "/", mustWork = FALSE)
}

script_path <- detect_script_path()
script_dir <- dirname(script_path)
repo_root <- if (basename(script_dir) == "01_simulation") dirname(script_dir) else getwd()
setwd(repo_root)

run_all <- file.path("01_simulation", "run_all.R")
if (!file.exists(run_all)) {
  stop("Could not find ", run_all, ". Run this script from the repository root or from 01_simulation/.")
}

dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

steps <- data.frame(
  step = c(
    "main_pilot",
    "main_validate",
    "main_figures",
    "multicov_pilot",
    "multicov_validate",
    "multicov_figures",
    "nonlinear_pilot",
    "nonlinear_validate",
    "nonlinear_figures",
    "trunc_diagnose",
    "trunc_summarize",
    "trunc_figures"
  ),
  mode = c(
    "pilot-main",
    "validate-main",
    "make-figures",
    "pilot-multicov",
    "validate-multicov",
    "make-multicov-figures",
    "pilot-nonlinear",
    "validate-nonlinear",
    "make-nonlinear-figures",
    "diagnose-ripd-truncation",
    "summarize-ripd-truncation",
    "make-ripd-truncation-figures"
  ),
  output_root = c(
    rep(OUTPUT_ROOT_MAIN, 3),
    rep(OUTPUT_ROOT_MULTICOV, 3),
    rep(OUTPUT_ROOT_NONLINEAR, 3),
    rep(OUTPUT_ROOT_TRUNC, 3)
  ),
  use_workers = c(
    TRUE, FALSE, FALSE,
    TRUE, FALSE, FALSE,
    TRUE, FALSE, FALSE,
    TRUE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

tail_lines <- function(path, n = 50L) {
  if (!file.exists(path)) return(character())
  lines <- readLines(path, warn = FALSE)
  utils::tail(lines, n)
}

run_step <- function(step_row) {
  log_file <- file.path(LOG_DIR, paste0(step_row$step, ".log"))
  args <- c(
    run_all,
    paste0("--mode=", step_row$mode),
    paste0("--nsim=", NSIM),
    paste0("--output-root=", step_row$output_root)
  )
  n_workers_for_step <- NA_integer_
  if (isTRUE(step_row$use_workers)) {
    args <- c(args, paste0("--n-workers=", N_WORKERS))
    n_workers_for_step <- N_WORKERS
  }

  started_at <- Sys.time()
  message(sprintf("[%s] starting mode=%s output_root=%s log=%s",
                  step_row$step, step_row$mode, step_row$output_root, log_file))
  status <- system2(rscript, args = args, stdout = log_file, stderr = log_file)
  finished_at <- Sys.time()
  elapsed <- as.numeric(difftime(finished_at, started_at, units = "secs"))
  if (is.null(status)) status <- 0L
  status <- as.integer(status)

  message(sprintf("[%s] exit_status=%d elapsed_seconds=%.1f log=%s",
                  step_row$step, status, elapsed, log_file))

  result <- data.frame(
    step = step_row$step,
    mode = step_row$mode,
    nsim = NSIM,
    n_workers = n_workers_for_step,
    output_root = step_row$output_root,
    log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
    exit_status = status,
    elapsed_seconds = elapsed,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )

  if (status != 0L) {
    message("Failed command:")
    message(paste(c(shQuote(rscript), shQuote(args)), collapse = " "))
    message("Log file: ", normalizePath(log_file, winslash = "/", mustWork = FALSE))
    message("Last 50 log lines:")
    message(paste(tail_lines(log_file, 50L), collapse = "\n"))
    attr(result, "failed") <- TRUE
  }

  result
}

manifest_rows <- list()
failed <- FALSE
for (i in seq_len(nrow(steps))) {
  row <- run_step(steps[i, , drop = FALSE])
  manifest_rows[[length(manifest_rows) + 1L]] <- row
  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(LOG_DIR, "final_simulation_manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
  if (isTRUE(attr(row, "failed"))) {
    failed <- TRUE
    break
  }
}

manifest_path <- normalizePath(file.path(LOG_DIR, "final_simulation_manifest.csv"), winslash = "/", mustWork = FALSE)
if (failed) {
  stop("Final simulation workflow stopped after a failed step. Manifest: ", manifest_path)
}

message("Final simulation workflow completed successfully.")
message("Manifest: ", manifest_path)
