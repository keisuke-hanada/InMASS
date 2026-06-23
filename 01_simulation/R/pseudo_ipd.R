reconstruct_pseudo_ipd <- function(data_mean, data_var, formula_ma, meta_fit, strata, seed) {
  set.seed(seed)
  if (!isTRUE(meta_fit$converged)) return(data.frame())

  formula_ma <- stats::as.formula(formula_ma)
  vars <- all.vars(formula_ma)
  arm_var <- vars[2]
  covariates <- vars[-c(1, 2)]
  coefficients <- meta_fit$coefficients
  out <- vector("list", strata * 2L)
  pos <- 1L

  for (st in seq_len(strata)) {
    ad_mean <- data_mean[data_mean$strata == st, , drop = FALSE]
    ad_var <- data_var[data_var$strata == st, , drop = FALSE]
    arms <- sort(unique(ad_mean[[arm_var]]))

    for (arm in arms) {
      mean_row <- ad_mean[ad_mean[[arm_var]] == arm, , drop = FALSE]
      var_row <- ad_var[ad_var[[arm_var]] == arm, , drop = FALSE]
      n_arm <- as.integer(mean_row$n[1])
      ripd <- mean_row[rep(1, n_arm), , drop = FALSE]

      for (covariate in covariates) {
        ripd[[covariate]] <- stats::rnorm(
          n_arm,
          mean = mean_row[[covariate]][1],
          sd = sqrt(var_row[[covariate]][1])
        )
      }

      model_data <- stats::model.frame(formula_ma, data = ripd)
      y_hat <- as.numeric(stats::model.matrix(formula_ma, data = model_data) %*% coefficients)
      sd_y <- sqrt(max(var_row[[vars[1]]][1] - stats::var(y_hat), 0))
      model_data[[vars[1]]] <- y_hat + stats::rnorm(n_arm, sd = sd_y)
      model_data$strata <- st
      model_data$vi <- sd_y^2
      out[[pos]] <- model_data
      pos <- pos + 1L
    }
  }

  out <- out[seq_len(pos - 1L)]
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}
