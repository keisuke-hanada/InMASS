make_aggregate_data <- function(ipd, study_col = "strata", arm_col = "x1k") {
  split_key <- interaction(ipd[[study_col]], ipd[[arm_col]], drop = TRUE)
  groups <- split(ipd, split_key)

  rows <- lapply(groups, function(d) {
    numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    means <- vapply(d[numeric_cols], mean, numeric(1))
    vars <- vapply(d[numeric_cols], stats::var, numeric(1))
    means[study_col] <- d[[study_col]][1]
    vars[study_col] <- d[[study_col]][1]
    means[arm_col] <- d[[arm_col]][1]
    vars[arm_col] <- d[[arm_col]][1]
    list(
      mean = data.frame(as.list(means), var = "mean", n = nrow(d), check.names = FALSE),
      var = data.frame(as.list(vars), var = "var", n = nrow(d), check.names = FALSE)
    )
  })

  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  row.names(out) <- NULL
  out
}

target_aggregate_rows <- function(target_ipd, template_names, replicate, arm_col = "x1k") {
  target <- target_ipd
  target$strata <- 0
  target$nsim <- replicate
  ad <- make_aggregate_data(target, arm_col = arm_col)
  missing <- setdiff(template_names, names(ad))
  for (col in missing) ad[[col]] <- NA
  ad[, template_names, drop = FALSE]
}
