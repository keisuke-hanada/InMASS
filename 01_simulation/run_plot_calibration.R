# After running run_all.R, run: Rscript 01_simulation/run_plot_calibration.R
# This redraws only the nonlinear robustness MSE plots; it does not rerun simulations.

OUTPUT_ROOT <- NULL
# If OUTPUT_ROOT is NULL, the script uses the most recently modified
# robustness_nonlinear_results_nsim*.csv file found under the repository root.
# To force a specific run, set for example:
# OUTPUT_ROOT <- "results_final_nonlinear"

detect_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[startsWith(args, "--file=")]
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
    script_dir <- dirname(script_path)
    if (basename(script_dir) == "01_simulation") {
      return(dirname(script_dir))
    }
    return(script_dir)
  }
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (basename(wd) == "01_simulation") {
    return(dirname(wd))
  }
  wd
}

repo_root <- detect_repo_root()
setwd(repo_root)

plotting_file <- file.path("01_simulation", "R", "plotting.R")
if (!file.exists(plotting_file)) {
  stop("Could not find ", plotting_file, ". Run from the repository root.")
}
source(plotting_file)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for figure generation.")
}

find_nonlinear_results <- function(output_root = NULL) {
  pattern <- "^robustness_nonlinear_results_nsim[0-9]+\\.csv$"
  if (!is.null(output_root) && nzchar(output_root)) {
    candidates <- list.files(file.path(output_root, "summary"), pattern = pattern, full.names = TRUE)
  } else {
    candidates <- list.files(".", pattern = pattern, recursive = TRUE, full.names = TRUE)
  }
  if (!length(candidates)) {
    stop("No robustness_nonlinear_results_nsim*.csv file found.")
  }
  info <- file.info(candidates)
  candidates[order(info$mtime, decreasing = TRUE)][1]
}

results_file <- find_nonlinear_results(OUTPUT_ROOT)
output_root <- dirname(dirname(results_file))
out_root <- file.path(output_root, "figures", "supplement", "robustness_nonlinear")
if (!dir.exists(out_root)) {
  stop("Expected output directory does not exist: ", out_root)
}

results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
d <- prepare_legacy_plot_data(results)
d$model <- "nonlinear"
d$Method <- d$Method_main
d$formula <- factor(d$formula, levels = c(2, 1))

labs <- legacy_labellers(d)
labs$model <- c("Nonlinear" = "nonlinear")

mse_data <- filter_mse_power_methods(aggregate_mse_legacy(d, "Method_main"))
mse_data$formula <- factor(mse_data$formula, levels = c(2, 1))

mse_without_3to1 <- mse_data[as.character(mse_data$Method) != "3:1", , drop = FALSE]
mse_without_3to1$Method <- droplevels(mse_without_3to1$Method)

mse_3to1_only <- mse_data[as.character(mse_data$Method) == "3:1", , drop = FALSE]
mse_3to1_only$Method <- droplevels(mse_3to1_only$Method)

if (!nrow(mse_without_3to1)) {
  stop("No nonlinear MSE rows remain after excluding Method == '3:1'.")
}
if (!nrow(mse_3to1_only)) {
  stop("No nonlinear MSE rows found for Method == '3:1'.")
}
if (any(as.character(mse_without_3to1$Method) == "3:1")) {
  stop("Filtering failed: Method == '3:1' remains in the main nonlinear MSE plot data.")
}
if (!all(as.character(mse_3to1_only$Method) == "3:1")) {
  stop("Filtering failed: the 3:1-only plot data contain other methods.")
}

main_path <- file.path(out_root, "robustness_nonlinear_mse.pdf")
only_path <- file.path(out_root, "robustness_nonlinear_mse_3to1_only.pdf")

save_plot(
  legacy_mse_plot(mse_without_3to1, labs, stats::as.formula("formula ~ strata")),
  main_path,
  7,
  4.8
)

plot_3to1_only <- legacy_mse_plot(mse_3to1_only, labs, stats::as.formula("formula ~ strata")) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

save_plot(
  plot_3to1_only,
  only_path,
  7,
  4.8
)

message("Input results: ", normalizePath(results_file, winslash = "/", mustWork = TRUE))
message("Wrote: ", normalizePath(main_path, winslash = "/", mustWork = TRUE))
message("Wrote: ", normalizePath(only_path, winslash = "/", mustWork = TRUE))
