load_v1_summary <- function(paths) {
  subfolders <- c("test_1to1" = "1to1", "test_3to1" = "3to1", "test_4to0" = "4to0")
  files <- file.path(paths$v1_root, names(subfolders), "total_evaluates.csv")
  existing <- file.exists(files)
  if (!any(existing)) return(data.frame())

  rows <- lapply(which(existing), function(i) {
    d <- utils::read.csv(files[i], stringsAsFactors = FALSE)
    d$allocation <- unname(subfolders[[basename(dirname(files[i]))]])
    d
  })
  do.call(rbind, rows)
}

make_v1_comparison_table <- function(v1) {
  if (!nrow(v1)) return(data.frame())
  v1$formula_id <- ifelse(v1$formula == 1, "misspecified", "correct")
  v1$estimator <- ifelse(
    v1$method == "ipd_method", "target_only",
    ifelse(v1$method == "propose_method", "inmass", "meta_regression")
  )
  v1$covariate_distribution <- ifelse(v1$model == "make_model1", "normal", "chi2")
  v1$K <- v1$strata

  mean_rows <- v1[v1$eval == "meanfunc", ]
  var_rows <- v1[v1$eval == "varfunc", ]
  bias_rows <- v1[v1$eval == "biasfunc", ]

  merge_keys <- c("allocation", "K", "n", "covariate_distribution", "estimator", "formula_id")
  mse_rows <- merge(
    mean_rows[c(merge_keys, "mean")],
    var_rows[c(merge_keys, "mean")],
    by = merge_keys,
    suffixes = c("_estimate", "_variance")
  )
  mse_rows$v1_mse <- mse_rows$mean_variance + (mse_rows$mean_estimate - 2)^2

  bias_rows <- bias_rows[c(merge_keys, "mean")]
  names(bias_rows)[names(bias_rows) == "mean"] <- "v1_bias"
  merge(mse_rows[c(merge_keys, "v1_mse")], bias_rows, by = merge_keys, all = TRUE)
}

validate_pilot_against_v1 <- function(paths, nsim = 10L) {
  v2_file <- file.path(paths$summary, sprintf("main_summary_nsim%d.csv", nsim))
  if (!file.exists(v2_file)) {
    stop("Run the v2 pilot before validation: missing ", v2_file)
  }
  v2 <- utils::read.csv(v2_file, stringsAsFactors = FALSE)
  v1 <- make_v1_comparison_table(load_v1_summary(paths))
  if (!nrow(v1)) {
    message("No v1 generated summaries found under ", paths$v1_root, "; writing v2-only validation summary.")
    out <- v2
    out$validation_note <- "v1 outputs not found"
  } else {
    keys <- c("allocation", "K", "n", "covariate_distribution", "estimator", "formula_id")
    out <- merge(v2, v1, by = keys, all.x = TRUE)
    out$bias_difference_v2_minus_v1 <- out$bias - out$v1_bias
    out$mse_difference_v2_minus_v1 <- out$mse - out$v1_mse
    out$validation_note <- ifelse(is.na(out$v1_bias) & is.na(out$v1_mse), "no matching v1 row", "matched")
  }
  write_csv(out, file.path(paths$summary, sprintf("pilot_validation_against_v1_nsim%d.csv", nsim)))
  out
}
