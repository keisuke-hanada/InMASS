generate_covariate <- function(n, distribution, mean_shift = 0) {
  if (distribution == "normal") {
    stats::rnorm(n) + mean_shift
  } else if (distribution == "chi2") {
    ((stats::rnorm(n)^2 + stats::rnorm(n)^2) / 2) + mean_shift - 1
  } else {
    stop("Unknown covariate distribution: ", distribution)
  }
}

target_allocation <- function(n, allocation) {
  if (allocation == "1to1") {
    rep(c(0, 1), length.out = n)
  } else if (allocation == "3to1") {
    rep(c(0, 1, 1, 1), length.out = n)
  } else if (allocation == "4to0") {
    rep(1, n)
  } else {
    stop("Unknown allocation: ", allocation)
  }
}

external_allocation <- function(n) {
  c(rep(0, n / 2), rep(1, n / 2))
}

simulate_outcome_main <- function(z, x, sigma) {
  beta0 <- 1
  delta_t <- 2
  beta_x <- -1
  beta_zx <- 0.5
  beta0 + delta_t * z + beta_x * x + beta_zx * z * x + stats::rnorm(length(z), sd = sigma)
}

simulate_outcome_multicov <- function(z, x1, x2, sigma) {
  beta0 <- 1
  delta_t <- 2
  beta1 <- -1
  beta2 <- 0.5
  beta3 <- 0.5
  beta4 <- -0.25
  beta0 + delta_t * z + beta1 * x1 + beta2 * x2 + beta3 * z * x1 + beta4 * z * x2 +
    stats::rnorm(length(z), sd = sigma)
}

generate_main_data <- function(spec, seed) {
  set.seed(seed)
  K <- as.integer(spec$K)
  n <- as.integer(spec$n)
  nsim <- as.integer(spec$nsim)
  sigma <- spec$sigma
  distribution <- spec$covariate_distribution

  n_external <- 2L * round(stats::runif(K, n / 2, 2 * n))
  mu_external <- if (K == 1L) 0 else 4 * (seq_len(K) - 1) / (K - 1) - 1

  target_list <- vector("list", nsim)
  external_list <- vector("list", nsim)
  aggregate_list <- vector("list", nsim)

  for (replicate in seq_len(nsim)) {
    z_t <- target_allocation(n, spec$allocation)
    x_t <- generate_covariate(n, distribution, 0)
    target_ipd <- data.frame(
      x1k = z_t,
      x2k = x_t,
      x3k = generate_covariate(n, "normal", 0),
      yik = simulate_outcome_main(z_t, x_t, sigma),
      nsim = replicate
    )

    external_by_study <- lapply(seq_len(K), function(k) {
      nk <- n_external[k]
      z_k <- external_allocation(nk)
      x_k <- generate_covariate(nk, distribution, mu_external[k])
      data.frame(
        x1k = z_k,
        x2k = x_k,
        x3k = generate_covariate(nk, "normal", 0),
        yik = simulate_outcome_main(z_k, x_k, sigma),
        strata = k,
        nsim = replicate
      )
    })

    external_ipd <- do.call(rbind, external_by_study)
    row.names(external_ipd) <- NULL
    aggregate_data <- make_aggregate_data(external_ipd)
    aggregate_data$nsim <- replicate

    target_list[[replicate]] <- target_ipd
    external_list[[replicate]] <- external_ipd
    aggregate_list[[replicate]] <- aggregate_data
  }

  params <- list(
    scenario_id = spec$scenario_id,
    allocation = spec$allocation,
    n = n,
    K = K,
    strata = K,
    nsim = nsim,
    sigma = sigma,
    dist_x2 = distribution,
    formula = spec$formula_ma,
    treatment_var = "x1k",
    simulation_family = spec$simulation_family %||% "main",
    dgm = spec$dgm,
    beta = c(1, 2, -1, 0, 0.5),
    truth = spec$truth,
    n_external = n_external
  )

  list(
    target_ipd = do.call(rbind, target_list),
    strata_ipd = do.call(rbind, external_list),
    strata_ad = do.call(rbind, aggregate_list),
    params = params
  )
}

generate_multicov_data <- function(spec, seed) {
  set.seed(seed)
  K <- as.integer(spec$K)
  n <- as.integer(spec$n)
  nsim <- as.integer(spec$nsim)
  sigma <- spec$sigma

  n_external <- 2L * round(stats::runif(K, n / 2, 2 * n))
  mu1_external <- if (K == 1L) 0 else 4 * (seq_len(K) - 1) / (K - 1) - 1
  mu2_external <- if (K == 1L) 0 else 2 * (seq_len(K) - 1) / (K - 1) - 0.5

  target_list <- vector("list", nsim)
  external_list <- vector("list", nsim)
  aggregate_list <- vector("list", nsim)

  for (replicate in seq_len(nsim)) {
    z_t <- target_allocation(n, spec$allocation)
    x1_t <- stats::rnorm(n, mean = 0, sd = 1)
    x2_t <- stats::rnorm(n, mean = 0, sd = 1)
    target_ipd <- data.frame(
      z = z_t,
      x1 = x1_t,
      x2 = x2_t,
      yik = simulate_outcome_multicov(z_t, x1_t, x2_t, sigma),
      nsim = replicate
    )

    external_by_study <- lapply(seq_len(K), function(k) {
      nk <- n_external[k]
      z_k <- external_allocation(nk)
      x1_k <- stats::rnorm(nk, mean = mu1_external[k], sd = 1)
      x2_k <- stats::rnorm(nk, mean = mu2_external[k], sd = 1)
      data.frame(
        z = z_k,
        x1 = x1_k,
        x2 = x2_k,
        yik = simulate_outcome_multicov(z_k, x1_k, x2_k, sigma),
        strata = k,
        nsim = replicate
      )
    })

    external_ipd <- do.call(rbind, external_by_study)
    row.names(external_ipd) <- NULL
    aggregate_data <- make_aggregate_data(external_ipd, arm_col = "z")
    aggregate_data$nsim <- replicate

    target_list[[replicate]] <- target_ipd
    external_list[[replicate]] <- external_ipd
    aggregate_list[[replicate]] <- aggregate_data
  }

  params <- list(
    scenario_id = spec$scenario_id,
    allocation = spec$allocation,
    n = n,
    K = K,
    strata = K,
    nsim = nsim,
    sigma = sigma,
    formula = spec$formula_ma,
    treatment_var = "z",
    simulation_family = "robustness_multicov",
    dgm = "multicov",
    beta = c(delta_T = 2, beta0 = 1, beta1 = -1, beta2 = 0.5, beta3 = 0.5, beta4 = -0.25),
    truth = spec$truth,
    n_external = n_external,
    mu1_external = mu1_external,
    mu2_external = mu2_external
  )

  list(
    target_ipd = do.call(rbind, target_list),
    strata_ipd = do.call(rbind, external_list),
    strata_ad = do.call(rbind, aggregate_list),
    params = params
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
