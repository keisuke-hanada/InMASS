evaluate_estimates <- function(results, truth = 2, null_value = 0) {
  keys <- unique(results[c("scenario_id", "estimator", "formula_id")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    idx <- results$scenario_id == key$scenario_id &
      results$estimator == key$estimator &
      results$formula_id == key$formula_id
    d <- results[idx & results$converged, , drop = FALSE]

    if (!nrow(d)) {
      return(data.frame(
        key,
        n_replicates = sum(idx),
        n_converged = 0L,
        bias = NA_real_,
        mse = NA_real_,
        power = NA_real_,
        mean_estimate = NA_real_,
        mean_se = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      key,
      n_replicates = sum(idx),
      n_converged = nrow(d),
      bias = mean(d$estimate - truth, na.rm = TRUE),
      mse = mean((d$estimate - truth)^2, na.rm = TRUE),
      power = mean(d$ci_low > null_value | d$ci_high < null_value, na.rm = TRUE),
      mean_estimate = mean(d$estimate, na.rm = TRUE),
      mean_se = mean(d$se, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
