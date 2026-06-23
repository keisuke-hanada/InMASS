estimate_density_ratio <- function(target_ipd, pseudo_ipd, formula_ma, ps_meta = 1) {
  formula_ma <- stats::as.formula(formula_ma)
  vars <- all.vars(formula_ma)
  arm_var <- vars[2]

  target_ipd$vi <- NA_real_
  if (any(target_ipd[[arm_var]] == 0)) {
    target_ipd$vi[target_ipd[[arm_var]] == 0] <- stats::var(target_ipd[[vars[1]]][target_ipd[[arm_var]] == 0])
  }
  if (any(target_ipd[[arm_var]] == 1)) {
    target_ipd$vi[target_ipd[[arm_var]] == 1] <- stats::var(target_ipd[[vars[1]]][target_ipd[[arm_var]] == 1])
  }
  target_ipd$strata <- 0
  target_aligned <- target_ipd
  missing <- setdiff(names(pseudo_ipd), names(target_aligned))
  for (col in missing) target_aligned[[col]] <- NA
  target_aligned <- target_aligned[, names(pseudo_ipd), drop = FALSE]

  if (ps_meta == 1) {
    combined <- rbind(target_aligned, pseudo_ipd)
    source_count_mode <- "all_rows_legacy"
  } else if (ps_meta == 2) {
    combined <- rbind(target_aligned, pseudo_ipd[pseudo_ipd[[arm_var]] == 0, , drop = FALSE])
    source_count_mode <- "source_rows"
  } else {
    stop("ps_meta must be 1 or 2.")
  }

  combined$set <- as.numeric(combined$strata == 0)
  covariates <- vars[-c(1, 2)]
  quadratic <- paste0("I(", covariates, "^2)")
  ps_formula <- stats::reformulate(c(covariates, quadratic), response = "id")
  ps_dat <- rbind(data.frame(id = 1, combined[combined$set == 1, , drop = FALSE]),
                  data.frame(id = 0, combined))
  ps_model <- stats::glm(ps_formula, family = stats::binomial(link = "logit"), data = ps_dat)
  psval <- stats::predict(ps_model, type = "response")[ps_dat$id == 0]
  n_target <- sum(ps_dat$id == 1)
  n_source <- if (ps_meta == 1) length(ps_dat$id == 0) else sum(ps_dat$id == 0)
  combined$weight <- psval / (1 - psval) * (n_source / n_target)
  if (ps_meta == 2) combined$weight[combined[[arm_var]] == 1] <- 1

  list(
    data = combined,
    ps_model = ps_model,
    diagnostics = list(
      ps_meta = ps_meta,
      source_count_mode = source_count_mode,
      density_ratio_covariates = paste(covariates, collapse = "+"),
      density_ratio_quadratic_terms = paste(quadratic, collapse = "+")
    )
  )
}
