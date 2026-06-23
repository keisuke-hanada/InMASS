reconstruct_pseudo_ipd <- function(data_mean, data_var, formula_ma, meta_fit, strata, seed,
                                   metadata = list()) {
  set.seed(seed)
  if (!isTRUE(meta_fit$converged)) {
    out <- data.frame()
    attr(out, "truncation_diagnostics") <- data.frame()
    return(out)
  }

  formula_ma <- stats::as.formula(formula_ma)
  vars <- all.vars(formula_ma)
  arm_var <- vars[2]
  covariates <- vars[-c(1, 2)]
  derived_covariates <- intersect(c("x_second", "x_centered_second"), covariates)
  base_covariates <- setdiff(covariates, derived_covariates)
  coefficients <- meta_fit$coefficients
  out <- vector("list", strata * 2L)
  truncation <- vector("list", strata * 2L)
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

      for (covariate in base_covariates) {
        ripd[[covariate]] <- stats::rnorm(
          n_arm,
          mean = mean_row[[covariate]][1],
          sd = sqrt(var_row[[covariate]][1])
        )
      }
      if ("x_second" %in% derived_covariates && "x" %in% names(ripd)) {
        ripd$x_second <- ripd$x^2
      }
      if ("x_centered_second" %in% derived_covariates && "x" %in% names(ripd)) {
        ripd$x_centered_second <- ripd$x^2 - 1
      }

      model_data <- stats::model.frame(formula_ma, data = ripd)
      y_hat <- as.numeric(stats::model.matrix(formula_ma, data = model_data) %*% coefficients)
      reported_outcome_variance <- var_row[[vars[1]]][1]
      model_explained_covariate_variance <- stats::var(y_hat)
      resid_var_untruncated <- reported_outcome_variance - model_explained_covariate_variance
      resid_var_used <- max(resid_var_untruncated, 0)
      sd_y <- sqrt(resid_var_used)
      model_data[[vars[1]]] <- y_hat + stats::rnorm(n_arm, sd = sd_y)
      model_data$strata <- st
      model_data$vi <- sd_y^2
      out[[pos]] <- model_data
      truncation[[pos]] <- c(
        metadata,
        list(
          study_id = st,
          arm = arm,
          reported_outcome_variance = reported_outcome_variance,
          model_explained_covariate_variance = model_explained_covariate_variance,
          resid_var_untruncated = resid_var_untruncated,
          resid_var_used = resid_var_used,
          truncated = resid_var_untruncated <= 0
        )
      )
      pos <- pos + 1L
    }
  }

  out <- out[seq_len(pos - 1L)]
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  truncation <- truncation[seq_len(pos - 1L)]
  truncation <- do.call(rbind, lapply(truncation, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  row.names(truncation) <- NULL
  numeric_cols <- c(
    "replicate", "K", "n", "study_id", "arm", "reported_outcome_variance",
    "model_explained_covariate_variance", "resid_var_untruncated", "resid_var_used"
  )
  for (col in intersect(numeric_cols, names(truncation))) {
    truncation[[col]] <- as.numeric(truncation[[col]])
  }
  if ("truncated" %in% names(truncation)) {
    truncation$truncated <- as.logical(truncation$truncated)
  }
  attr(out, "truncation_diagnostics") <- truncation
  out
}
