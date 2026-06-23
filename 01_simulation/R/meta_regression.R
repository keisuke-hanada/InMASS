fit_meta_regression <- function(data_mean, data_var, formula, method = "DL") {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for meta-regression.")
  }

  formula <- stats::as.formula(formula)
  yi <- data_mean$yik
  vi <- data_var$yik
  fit <- tryCatch(
    metafor::rma.uni(
      yi = yi,
      vi = vi,
      mods = stats::update(formula, NULL ~ .),
      data = data_mean,
      method = method
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(converged = FALSE, error = conditionMessage(fit)))
  }

  list(
    converged = TRUE,
    fit = fit,
    coefficients = as.numeric(fit$b),
    coefficient_names = row.names(fit$b),
    vcov = stats::vcov(fit),
    se = as.numeric(fit$se),
    ci_low = as.numeric(fit$ci.lb),
    ci_high = as.numeric(fit$ci.ub),
    tau2 = if (!is.null(fit$tau2)) as.numeric(fit$tau2) else NA_real_,
    k = if (!is.null(fit$k)) as.integer(fit$k) else NA_integer_
  )
}

normalize_meta_coefficient_names <- function(names) {
  names[names == "(Intercept)"] <- "intrcpt"
  names
}

extract_meta_coefficient <- function(meta_fit, coefficient = "x1k") {
  if (!isTRUE(meta_fit$converged)) {
    return(list(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  idx <- match(coefficient, meta_fit$coefficient_names)
  if (is.na(idx)) {
    return(list(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  list(
    estimate = meta_fit$coefficients[idx],
    se = meta_fit$se[idx],
    ci_low = meta_fit$ci_low[idx],
    ci_high = meta_fit$ci_high[idx]
  )
}

plugin_contrast <- function(meta_fit, formula, target_ipd, level = 0.95) {
  if (!isTRUE(meta_fit$converged)) {
    return(list(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      converged = FALSE, diagnostics = meta_fit$error %||% "meta-regression failed"
    ))
  }

  formula <- stats::as.formula(formula)
  vars <- all.vars(formula)
  response <- vars[1]
  arm_var <- vars[2]
  covariates <- setdiff(vars[-c(1, 2)], response)

  arm_profiles <- lapply(c(0, 1), function(arm) {
    d <- target_ipd[target_ipd[[arm_var]] == arm, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    profile <- target_ipd[1, , drop = FALSE]
    profile[] <- NA
    profile[[response]] <- 0
    profile[[arm_var]] <- arm
    for (covariate in covariates) {
      profile[[covariate]] <- mean(d[[covariate]], na.rm = TRUE)
    }
    profile
  })
  names(arm_profiles) <- c("control", "treated")

  profile_note <- "target arm-specific covariate means used"
  if (is.null(arm_profiles$treated)) {
    return(list(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      converged = FALSE, diagnostics = "treated target profile unavailable"
    ))
  }
  if (is.null(arm_profiles$control)) {
    return(list(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      converged = FALSE,
      diagnostics = "plugin estimator not estimable: target control covariate summaries are unavailable"
    ))
  }

  profiles <- rbind(arm_profiles$control, arm_profiles$treated)
  mm <- stats::model.matrix(formula, data = stats::model.frame(formula, profiles, na.action = stats::na.pass))
  contrast_model <- mm[2, ] - mm[1, ]
  names(contrast_model) <- normalize_meta_coefficient_names(names(contrast_model))

  coefficient_names <- meta_fit$coefficient_names
  contrast <- stats::setNames(rep(0, length(coefficient_names)), coefficient_names)
  common <- intersect(names(contrast_model), coefficient_names)
  contrast[common] <- contrast_model[common]
  missing <- setdiff(names(contrast_model)[contrast_model != 0], coefficient_names)
  if (length(missing)) {
    return(list(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      converged = FALSE,
      diagnostics = paste("contrast coefficient names missing from meta-regression fit:", paste(missing, collapse = ","))
    ))
  }

  beta <- stats::setNames(meta_fit$coefficients, coefficient_names)
  vcov <- meta_fit$vcov
  if (is.null(dimnames(vcov))) dimnames(vcov) <- list(coefficient_names, coefficient_names)
  vcov <- vcov[coefficient_names, coefficient_names, drop = FALSE]
  estimate <- sum(contrast * beta)
  variance <- as.numeric(t(contrast) %*% vcov %*% contrast)
  if (!is.finite(variance) || variance < 0) {
    return(list(
      estimate = estimate, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      converged = FALSE, diagnostics = "delta-method variance is negative or non-finite"
    ))
  }
  se <- sqrt(variance)
  z <- stats::qnorm(1 - (1 - level) / 2)
  nonzero <- names(contrast)[contrast != 0]
  diagnostics <- paste0(
    "mapping=arm-level predicted treated-control contrast;",
    "beta0_hat=", ifelse(arm_var %in% nonzero, arm_var, "none"), ";",
    "beta_M_hat=", paste(setdiff(nonzero, arm_var), collapse = "+"), ";",
    "contrast=", paste(paste(names(contrast), signif(contrast, 6), sep = "="), collapse = ","), ";",
    profile_note
  )

  list(
    estimate = estimate,
    se = se,
    ci_low = estimate - z * se,
    ci_high = estimate + z * se,
    converged = TRUE,
    diagnostics = diagnostics
  )
}
