main_analysis_formulas <- function() {
  data.frame(
    formula_id = c("misspecified", "correct"),
    formula_label = c("Model misspecified", "Model specified"),
    formula = c("yik ~ 1 + x1k", "yik ~ 1 + x1k + x2k + x1k:x2k"),
    stringsAsFactors = FALSE
  )
}

build_main_scenario_grid <- function(nsim = 10000) {
  grid <- expand.grid(
    allocation = c("1to1", "3to1", "4to0"),
    K = c(5L, 10L, 30L),
    n = c(20L, 40L, 100L),
    covariate_distribution = c("normal", "chi2"),
    stringsAsFactors = FALSE
  )

  grid$scenario_id <- sprintf(
    "main_%s_K%02d_n%03d_%s",
    grid$allocation,
    grid$K,
    grid$n,
    grid$covariate_distribution
  )
  grid$nsim <- as.integer(nsim)
  grid$sigma <- 1
  grid$dgm <- "main_one_covariate"
  grid$truth <- 2
  grid$formula_ma <- "yik ~ 1 + x1k + x2k + x1k:x2k"
  grid$treatment_var <- "x1k"
  grid$simulation_family <- "main"
  grid
}

build_main_run_grid <- function(nsim = 10000) {
  merge(build_main_scenario_grid(nsim), main_analysis_formulas(), by = NULL)
}

multicov_analysis_formulas <- function() {
  data.frame(
    formula_id = c("misspecified", "correct"),
    formula_label = c("Model misspecified", "Model specified"),
    formula = c("yik ~ 1 + z", "yik ~ 1 + z + x1 + x2 + z:x1 + z:x2"),
    stringsAsFactors = FALSE
  )
}

build_multicov_scenario_grid <- function(nsim = 10000) {
  grid <- expand.grid(
    allocation = c("1to1", "3to1", "4to0"),
    K = c(5L, 10L, 30L),
    n = c(20L, 40L, 100L),
    stringsAsFactors = FALSE
  )
  grid$scenario_id <- sprintf(
    "robustness_multicov_%s_K%02d_n%03d",
    grid$allocation,
    grid$K,
    grid$n
  )
  grid$nsim <- as.integer(nsim)
  grid$sigma <- 1
  grid$dgm <- "multicov"
  grid$truth <- 2
  grid$formula_ma <- "yik ~ 1 + z + x1 + x2 + z:x1 + z:x2"
  grid$treatment_var <- "z"
  grid$simulation_family <- "robustness_multicov"
  grid$covariate_distribution <- "multicov_normal"
  grid
}

build_multicov_run_grid <- function(nsim = 10000) {
  merge(build_multicov_scenario_grid(nsim), multicov_analysis_formulas(), by = NULL)
}

nonlinear_analysis_formulas <- function() {
  data.frame(
    formula_id = c("misspecified", "correct"),
    formula_label = c("Model misspecified", "Model specified"),
    formula = c(
      "yik ~ 1 + z + x + z:x",
      "yik ~ 1 + z + x + x_second + z:x + z:x_centered_second"
    ),
    stringsAsFactors = FALSE
  )
}

build_nonlinear_scenario_grid <- function(nsim = 10000) {
  grid <- expand.grid(
    allocation = c("1to1", "3to1", "4to0"),
    K = c(5L, 10L, 30L),
    n = c(20L, 40L, 100L),
    stringsAsFactors = FALSE
  )
  grid$scenario_id <- sprintf(
    "robustness_nonlinear_%s_K%02d_n%03d",
    grid$allocation,
    grid$K,
    grid$n
  )
  grid$nsim <- as.integer(nsim)
  grid$sigma <- 1
  grid$dgm <- "nonlinear"
  grid$truth <- 2
  grid$formula_ma <- "yik ~ 1 + z + x + x_second + z:x + z:x_centered_second"
  grid$treatment_var <- "z"
  grid$simulation_family <- "robustness_nonlinear"
  grid$covariate_distribution <- "nonlinear_normal"
  grid$density_covariates <- "x+x_second"
  grid$density_include_quadratic <- FALSE
  grid
}

build_nonlinear_run_grid <- function(nsim = 10000) {
  merge(build_nonlinear_scenario_grid(nsim), nonlinear_analysis_formulas(), by = NULL)
}
