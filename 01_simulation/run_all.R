source(file.path("01_simulation", "config", "paths.R"))
source(file.path("01_simulation", "config", "scenarios.R"))
source(file.path("01_simulation", "R", "seed_utils.R"))
source(file.path("01_simulation", "R", "io.R"))
source(file.path("01_simulation", "R", "aggregate_data.R"))
source(file.path("01_simulation", "R", "data_generation.R"))
source(file.path("01_simulation", "R", "meta_regression.R"))
source(file.path("01_simulation", "R", "pseudo_ipd.R"))
source(file.path("01_simulation", "R", "ripd_truncation.R"))
source(file.path("01_simulation", "R", "density_ratio.R"))
source(file.path("01_simulation", "R", "inmass.R"))
source(file.path("01_simulation", "R", "estimators.R"))
source(file.path("01_simulation", "R", "evaluation.R"))
source(file.path("01_simulation", "R", "plotting.R"))
source(file.path("01_simulation", "scripts", "run_main_scenarios.R"))
source(file.path("01_simulation", "scripts", "run_multicov_scenarios.R"))
source(file.path("01_simulation", "scripts", "run_nonlinear_scenarios.R"))
source(file.path("01_simulation", "scripts", "run_ripd_truncation_diagnostics.R"))
source(file.path("01_simulation", "scripts", "validate_pilot_internal.R"))
source(file.path("01_simulation", "scripts", "validate_multicov_internal.R"))
source(file.path("01_simulation", "scripts", "validate_nonlinear_internal.R"))
source(file.path("01_simulation", "scripts", "validate_pilot_against_v1.R"))
source(file.path("01_simulation", "scripts", "check_reproducibility.R"))
source(file.path("01_simulation", "scripts", "make_figures_tables.R"))

parse_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

mode <- parse_arg("mode", "pilot")
default_nsim <- if (mode %in% c("pilot", "validate-internal", "validate-pilot", "pilot-multicov", "validate-multicov", "pilot-nonlinear", "validate-nonlinear")) {
  "10"
} else if (mode %in% c("make-figures", "figures", "make-multicov-figures", "make-nonlinear-figures", "make-ripd-truncation-figures")) {
  "10"
} else if (mode %in% c("diagnose-ripd-truncation", "summarize-ripd-truncation")) {
  "10"
} else if (mode == "check-reproducibility") {
  "2"
} else {
  "10000"
}
nsim <- as.integer(parse_arg("nsim", default_nsim))
base_seed <- as.integer(parse_arg("base-seed", "1234"))
output_root <- parse_arg("output-root", "results")
paths <- simulation_paths(getwd(), output_root = output_root)

if (mode == "pilot") {
  invisible(run_main_scenarios(paths, nsim = nsim, base_seed = base_seed))
} else if (mode == "validate-internal") {
  invisible(validate_pilot_internal(paths, nsim = nsim))
} else if (mode == "validate-pilot") {
  invisible(validate_pilot_against_v1(paths, nsim = nsim))
} else if (mode == "pilot-multicov") {
  invisible(run_multicov_scenarios(paths, nsim = nsim, base_seed = base_seed))
} else if (mode == "validate-multicov") {
  invisible(validate_multicov_internal(paths, nsim = nsim))
} else if (mode == "pilot-nonlinear") {
  invisible(run_nonlinear_scenarios(paths, nsim = nsim, base_seed = base_seed))
} else if (mode == "validate-nonlinear") {
  invisible(validate_nonlinear_internal(paths, nsim = nsim))
} else if (mode == "diagnose-ripd-truncation") {
  invisible(run_ripd_truncation_diagnostics(paths, nsim = nsim, base_seed = base_seed))
} else if (mode == "summarize-ripd-truncation") {
  invisible(summarize_existing_ripd_truncation(paths, nsim = nsim))
} else if (mode == "check-reproducibility") {
  invisible(check_pilot_reproducibility(getwd(), nsim = nsim, base_seed = base_seed))
} else if (mode == "full") {
  invisible(run_main_scenarios(paths, nsim = nsim, base_seed = base_seed))
} else if (mode %in% c("make-figures", "figures")) {
  invisible(make_figures_tables(paths, nsim = nsim, pilot = output_root != "results"))
} else if (mode == "make-multicov-figures") {
  invisible(make_multicov_figures_tables(paths, nsim = nsim, pilot = output_root != "results"))
} else if (mode == "make-nonlinear-figures") {
  invisible(make_nonlinear_figures_tables(paths, nsim = nsim, pilot = output_root != "results"))
} else if (mode == "make-ripd-truncation-figures") {
  invisible(make_ripd_truncation_figures_tables(paths, nsim = nsim, pilot = output_root != "results"))
} else {
  stop("Unknown mode: ", mode)
}
