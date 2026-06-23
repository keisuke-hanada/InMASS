iwlm <- function(formula, data, weights) {
  y <- stats::model.response(stats::model.frame(formula, data))
  x <- stats::model.matrix(formula, data)
  w <- diag(weights)
  inv <- solve(t(x) %*% w %*% x)
  coef <- as.numeric(inv %*% t(x) %*% w %*% y)
  names(coef) <- colnames(x)
  residuals <- y - x %*% coef
  cov <- inv %*% (t(x) %*% w %*% diag(weights^2 * as.numeric(residuals)^2) %*% w %*% x) %*% inv
  dimnames(cov) <- list(colnames(x), colnames(x))
  out <- list(
    coefficients = coef,
    cov = cov,
    residuals = residuals,
    fitted.values = x %*% coef,
    formula = formula,
    data = data,
    weights = weights
  )
  class(out) <- "iwlm"
  out
}

coef.iwlm <- function(object, ...) {
  object$coefficients
}

confint.iwlm <- function(object, level = 0.95, ...) {
  x <- stats::model.matrix(object$formula, head(object$data))
  coef <- object$coefficients
  se <- sqrt(diag(object$cov))
  alpha <- 1 - level
  t_value <- stats::qt(1 - alpha / 2, df = sum(object$weights) - ncol(x))
  data.frame(Lower = coef - t_value * se, Upper = coef + t_value * se)
}

fit_inmass_core <- function(formula, formula_ma, data_mean, data_var, target_ipd, strata,
                            ps_meta, seed) {
  meta_fit <- fit_meta_regression(data_mean, data_var, formula_ma)
  pseudo <- reconstruct_pseudo_ipd(data_mean, data_var, formula_ma, meta_fit, strata, seed)
  if (!nrow(pseudo)) {
    return(list(converged = FALSE, error = "Pseudo-IPD reconstruction failed."))
  }
  weighted <- estimate_density_ratio(target_ipd, pseudo, formula_ma, ps_meta = ps_meta)
  fit <- tryCatch(
    iwlm(stats::as.formula(formula), data = weighted$data, weights = weighted$data$weight),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(converged = FALSE, error = conditionMessage(fit)))
  }
  list(converged = TRUE, fit = fit, meta_fit = meta_fit, diagnostics = weighted$diagnostics)
}
