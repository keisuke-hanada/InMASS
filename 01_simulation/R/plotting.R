figure_paths <- function(paths, pilot = TRUE) {
  if (pilot) {
    root <- file.path(dirname(paths$raw), "figures")
    list(main = file.path(root, "main"), supplement = file.path(root, "supplement"))
  } else {
    list(main = paths$figures_main, supplement = paths$figures_supplement)
  }
}

method_labels <- function(estimator) {
  labels <- c(
    target_only = "Target only",
    inmass = "InMASS",
    meta_regression = "Meta-regression",
    plugin = "Plug-in"
  )
  unname(labels[estimator])
}

legacy_model_name <- function(covariate_distribution) {
  ifelse(covariate_distribution == "normal", "make_model1", "make_model2")
}

legacy_formula_id <- function(formula_id) {
  ifelse(formula_id == "misspecified", 1L, 2L)
}

legacy_cformula <- function(formula) {
  ifelse(as.integer(as.character(formula)) == 1L, "MID", "ID")
}

legacy_method_suffix <- function(estimator) {
  suffix <- rep(NA_character_, length(estimator))
  suffix[estimator == "target_only"] <- ""
  suffix[estimator == "inmass"] <- "+InMASS"
  suffix[estimator == "meta_regression"] <- " MA"
  suffix[estimator == "plugin"] <- " Plug-in"
  suffix
}

legacy_allocation_prefix <- function(allocation, trailing_dash = TRUE) {
  prefix <- ifelse(allocation == "1to1", "1:1",
                   ifelse(allocation == "3to1", "3:1", "treatment"))
  if (trailing_dash) paste0(prefix, "-") else prefix
}

prepare_legacy_plot_data <- function(results) {
  d <- results
  d$estimate <- as.numeric(d$estimate)
  d$se <- as.numeric(d$se)
  d$ci_low <- as.numeric(d$ci_low)
  d$ci_high <- as.numeric(d$ci_high)
  d$true_delta <- as.numeric(d$true_delta)
  d$K <- as.integer(d$K)
  d$n <- as.integer(d$n)
  d$strata <- d$K
  d$model <- legacy_model_name(d$covariate_distribution)
  d$formula <- legacy_formula_id(d$formula_id)
  d$cformula <- legacy_cformula(d$formula)
  d$cmethod <- legacy_method_suffix(d$estimator)
  d$dname_dash <- legacy_allocation_prefix(d$allocation, trailing_dash = TRUE)
  d$dname <- legacy_allocation_prefix(d$allocation, trailing_dash = FALSE)
  d$Method_formula <- paste0(d$dname_dash, d$cformula, d$cmethod)
  d$Method_main <- paste0(d$dname, d$cmethod)
  d$bias <- d$estimate - d$true_delta
  d
}

aggregate_mse_legacy <- function(results, method_col = "Method_formula") {
  d <- results[results$converged & is.finite(results$estimate), , drop = FALSE]
  keys <- c("allocation", "strata", "n", "model", "formula", "cformula", method_col)
  groups <- split(d, interaction(d[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(g) {
    out <- g[1, keys, drop = FALSE]
    names(out)[names(out) == method_col] <- "Method"
    truth <- if ("true_delta" %in% names(g)) unique(g$true_delta) else 2
    truth <- truth[is.finite(truth)]
    if (length(truth) != 1L) {
      truth <- 2
    }
    out$mse <- mean((g$estimate - truth)^2, na.rm = TRUE)
    out
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

aggregate_power_curve_legacy <- function(results, method_col = "Method_formula", b.max = 4, h = 0.01) {
  d <- results[results$converged & is.finite(results$ci_low) & is.finite(results$ci_high), , drop = FALSE]
  keys <- c("allocation", "strata", "n", "model", "formula", "cformula", method_col)
  groups <- split(d, interaction(d[keys], drop = TRUE, lex.order = TRUE))
  b <- 2 + seq(0, (b.max - 2) / h) * h
  rows <- lapply(groups, function(g) {
    base <- g[1, keys, drop = FALSE]
    names(base)[names(base) == method_col] <- "Method"
    do.call(rbind, lapply(b, function(effect) {
      out <- base
      out$sd <- effect
      out$mean <- mean(effect <= g$ci_low | g$ci_high <= effect)
      out
    }))
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

validate_power_plot_data <- function(power_data, key_cols, path) {
  duplicate_keys <- duplicated(power_data[key_cols])
  has_replicate <- "replicate" %in% names(power_data)
  power_ok <- all(is.finite(power_data$mean) & power_data$mean >= 0 & power_data$mean <= 1)
  diagnostics <- data.frame(
    check = c("no_replicate_column", "no_duplicate_power_rows", "power_in_unit_interval"),
    passed = c(!has_replicate, !any(duplicate_keys), power_ok),
    failures = c(
      if (has_replicate) 1L else 0L,
      sum(duplicate_keys),
      sum(!(is.finite(power_data$mean) & power_data$mean >= 0 & power_data$mean <= 1))
    ),
    stringsAsFactors = FALSE
  )
  write_csv(diagnostics, path)
  if (!all(diagnostics$passed)) {
    stop("Power plotting data validation failed. See ", path)
  }
  diagnostics
}

bias_rows_legacy <- function(results, method_col = "Method_formula") {
  d <- results[results$converged & is.finite(results$bias), , drop = FALSE]
  out <- d[c("allocation", "strata", "n", "model", "formula", "cformula", method_col, "bias")]
  names(out)[names(out) == method_col] <- "Method"
  names(out)[names(out) == "bias"] <- "mean"
  row.names(out) <- NULL
  out
}

legacy_labellers <- function(d) {
  strata0 <- sort(unique(d$strata))
  n0 <- sort(unique(d$n))
  model0 <- c("make_model1", "make_model2")
  formula0 <- c(2, 1)
  strata.labs <- paste(strata0, " studies", sep = "")
  names(strata.labs) <- strata0
  n.labs <- paste("n=", n0, sep = "")
  names(n.labs) <- n0
  model.labs <- c("Normal", "Chi-squared")
  names(model.labs) <- model0
  formula.labs <- c("Model spesified", "Model misspesified")
  names(formula.labs) <- formula0
  list(strata = strata.labs, n = n.labs, model = model.labs, formula = formula.labs)
}

save_plot <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300)
  path
}

legacy_shape_scale <- function() {
  ggplot2::scale_shape_manual(values = rep(c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13), length.out = 32))
}

manuscript_mse_power_methods <- function() {
  c(
    "1:1",
    "1:1+InMASS",
    "1:1 Plug-in",
    "3:1",
    "3:1+InMASS",
    "3:1 Plug-in",
    "treatment+InMASS"
  )
}

supplement_1to1_methods <- function() {
  c(
    "1:1-ID",
    "1:1-ID+InMASS",
    "1:1-ID Plug-in",
    "1:1-MID",
    "1:1-MID+InMASS",
    "1:1-MID Plug-in"
  )
}

filter_mse_power_methods <- function(plot_data) {
  keep <- manuscript_mse_power_methods()
  plot_data <- plot_data[plot_data$Method %in% keep, , drop = FALSE]
  plot_data$Method <- factor(plot_data$Method, levels = keep)
  plot_data
}

filter_supplement_1to1_methods <- function(plot_data) {
  keep <- supplement_1to1_methods()
  plot_data <- plot_data[plot_data$Method %in% keep, , drop = FALSE]
  plot_data$Method <- factor(plot_data$Method, levels = keep)
  plot_data
}

validate_mse_power_plot_methods <- function(plot_data, path, plot_name) {
  allowed <- manuscript_mse_power_methods()
  methods <- unique(as.character(plot_data$Method))
  methods <- methods[!is.na(methods)]
  ma_methods <- grep(" MA$| MA |MA$", methods, value = TRUE)
  outside_allowed <- setdiff(methods, allowed)
  legend_levels <- levels(plot_data$Method)
  if (is.null(legend_levels)) legend_levels <- methods
  shown_missing_from_legend <- setdiff(methods, legend_levels)
  treatment_plugin <- grep("^treatment.*Plug-in", methods, value = TRUE)
  diagnostics <- data.frame(
    plot_name = plot_name,
    check = c(
      "no_meta_analysis_only_methods",
      "methods_subset_of_manuscript_set",
      "shown_methods_in_legend_levels",
      "no_treatment_only_plugin_method"
    ),
    passed = c(
      length(ma_methods) == 0L,
      length(outside_allowed) == 0L,
      length(shown_missing_from_legend) == 0L,
      length(treatment_plugin) == 0L
    ),
    failures = c(
      length(ma_methods),
      length(outside_allowed),
      length(shown_missing_from_legend),
      length(treatment_plugin)
    ),
    details = c(
      paste(ma_methods, collapse = "|"),
      paste(outside_allowed, collapse = "|"),
      paste(shown_missing_from_legend, collapse = "|"),
      paste(treatment_plugin, collapse = "|")
    ),
    stringsAsFactors = FALSE
  )
  write_csv(diagnostics, path)
  if (!all(diagnostics$passed)) {
    stop("MSE/Power plotting method validation failed for ", plot_name, ". See ", path)
  }
  diagnostics
}

append_plot_validation <- function(existing, plot_data, paths, nsim, plot_name, family = "main") {
  file <- file.path(paths$summary, sprintf("%s_mse_power_plot_method_validation_nsim%d.csv", family, nsim))
  current <- validate_mse_power_plot_methods(plot_data, tempfile(fileext = ".csv"), plot_name)
  if (is.null(existing)) {
    existing <- current
  } else {
    existing <- rbind(existing, current)
  }
  write_csv(existing, file)
  existing
}

validate_supplement_1to1_plot_methods <- function(plot_data, path, plot_name) {
  required <- supplement_1to1_methods()
  methods <- unique(as.character(plot_data$Method))
  methods <- methods[!is.na(methods)]
  old_labels <- intersect(methods, c("1:1", "1:1+InMASS", "1:1 Plug-in"))
  ma_methods <- grep(" MA$| MA |MA$", methods, value = TRUE)
  missing_required <- setdiff(required, methods)
  outside_required <- setdiff(methods, required)
  formula_levels <- sort(unique(as.character(plot_data$cformula)))
  legend_levels <- levels(plot_data$Method)
  if (is.null(legend_levels)) legend_levels <- methods
  legend_order_ok <- identical(legend_levels, required)
  diagnostics <- data.frame(
    plot_name = plot_name,
    check = c(
      "all_six_id_mid_labels_present",
      "no_old_undifferentiated_1to1_labels",
      "no_meta_analysis_only_methods",
      "both_model_specifications_included",
      "legend_levels_ordered"
    ),
    passed = c(
      length(missing_required) == 0L && length(outside_required) == 0L,
      length(old_labels) == 0L,
      length(ma_methods) == 0L,
      identical(formula_levels, c("ID", "MID")),
      legend_order_ok
    ),
    failures = c(
      length(missing_required) + length(outside_required),
      length(old_labels),
      length(ma_methods),
      length(setdiff(c("ID", "MID"), formula_levels)),
      if (legend_order_ok) 0L else 1L
    ),
    details = c(
      paste(c(paste0("missing:", missing_required), paste0("extra:", outside_required)), collapse = "|"),
      paste(old_labels, collapse = "|"),
      paste(ma_methods, collapse = "|"),
      paste(formula_levels, collapse = "|"),
      paste(legend_levels, collapse = "|")
    ),
    stringsAsFactors = FALSE
  )
  write_csv(diagnostics, path)
  if (!all(diagnostics$passed)) {
    stop("Supplement 1:1 plotting method validation failed for ", plot_name, ". See ", path)
  }
  diagnostics
}

append_supplement_1to1_validation <- function(existing, plot_data, paths, nsim, plot_name) {
  current <- validate_supplement_1to1_plot_methods(plot_data, tempfile(fileext = ".csv"), plot_name)
  if (is.null(existing)) {
    existing <- current
  } else {
    existing <- rbind(existing, current)
  }
  write_csv(existing, file.path(paths$summary, sprintf("main_supplement_1to1_plot_method_validation_nsim%d.csv", nsim)))
  existing
}

legacy_mse_plot <- function(g1dat, labs, facet = stats::as.formula("model ~ strata")) {
  ggplot2::ggplot(g1dat, ggplot2::aes(x = n, y = mse, col = Method, shape = Method, linetype = Method)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 0, intercept = 0) +
    legacy_shape_scale() +
    ggplot2::ylab("Mean Squared Error") +
    ggplot2::facet_grid(facet, labeller = ggplot2::labeller(strata = labs$strata, model = labs$model, formula = labs$formula)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = grid::unit(1.2, "cm"),
      plot.margin = ggplot2::margin(5.5, 5.5, 16, 5.5)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}

legacy_power_plot <- function(gdat, labs, title, facet = stats::as.formula("n ~ strata")) {
  ggplot2::ggplot(gdat, ggplot2::aes(x = sd - 2, y = mean, group = Method, col = Method, shape = Method, linetype = Method)) +
    ggplot2::geom_abline(slope = 0, intercept = 0.05, linetype = "dashed") +
    ggplot2::geom_abline(slope = 0, intercept = 0.8, linetype = "dashed") +
    ggplot2::geom_abline(slope = 0, intercept = 0.9, linetype = "dashed") +
    ggplot2::geom_line(linewidth = 1) +
    legacy_shape_scale() +
    ggplot2::ylab("Power") +
    ggplot2::xlab("Difference from CATE for target trial") +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::ggtitle(title) +
    ggplot2::facet_grid(facet, labeller = ggplot2::labeller(n = labs$n, strata = labs$strata, formula = labs$formula)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = grid::unit(1.2, "cm"),
      plot.margin = ggplot2::margin(5.5, 5.5, 16, 5.5)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}

legacy_bias_plot <- function(g4dat, labs, all_bias = FALSE, title = NULL) {
  p <- ggplot2::ggplot(g4dat, ggplot2::aes(x = Method, y = mean)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_boxplot() +
    ggplot2::ylab(expression("Bias of " ~ delta[T])) +
    ggplot2::facet_grid(formula ~ model, labeller = ggplot2::labeller(formula = labs$formula, model = labs$model)) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 12),
                   strip.text.x = ggplot2::element_text(size = 16),
                   strip.text.y = ggplot2::element_text(size = 16),
                   axis.title = ggplot2::element_text(size = 16)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1),
                    fill = ggplot2::guide_legend(nrow = 1))
  if (all_bias) {
    p <- p + ggplot2::coord_cartesian(ylim = c(-2, 2)) +
      ggplot2::labs(title = title, y = expression("Bias of " ~ delta[T]))
  }
  p
}

make_figures <- function(paths, nsim = 10L, pilot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for figure generation.")
  }
  results_file <- file.path(paths$summary, sprintf("main_results_nsim%d.csv", nsim))
  if (!file.exists(results_file)) stop("Replicate-level results file not found: ", results_file)
  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  d <- prepare_legacy_plot_data(results)
  labs <- legacy_labellers(d)
  out <- figure_paths(paths, pilot = pilot)
  dir.create(out$main, recursive = TRUE, showWarnings = FALSE)
  dir.create(out$supplement, recursive = TRUE, showWarnings = FALSE)
  generated <- character()
  method_validation <- NULL

  d_formula <- filter_supplement_1to1_methods(aggregate_mse_legacy(d, "Method_formula"))
  p_formula <- filter_supplement_1to1_methods(aggregate_power_curve_legacy(d, "Method_formula"))
  d_formula_control <- filter_mse_power_methods(aggregate_mse_legacy(d, "Method_main"))
  p_formula_control <- filter_mse_power_methods(aggregate_power_curve_legacy(d, "Method_main"))
  b_formula <- bias_rows_legacy(d, "Method_formula")

  d_1to1 <- d_formula[d_formula$allocation == "1to1", , drop = FALSE]
  p_1to1 <- p_formula[p_formula$allocation == "1to1", , drop = FALSE]
  supplement_1to1_validation <- NULL
  supplement_1to1_validation <- append_supplement_1to1_validation(supplement_1to1_validation, d_1to1, paths, nsim, "figure-1-1to1-mse")
  supplement_1to1_validation <- append_supplement_1to1_validation(supplement_1to1_validation, p_1to1, paths, nsim, "figure-1-1to1-power")
  generated <- c(generated,
    save_plot(legacy_mse_plot(d_1to1, labs), file.path(out$supplement, "figure-1-1to1-mse.pdf"), 7, 4.8),
    save_plot(legacy_power_plot(p_1to1[p_1to1$model == "make_model1", , drop = FALSE], labs, "Simulate model 1: normally distributed"),
              file.path(out$supplement, "figure-1-1to1-power-model1.pdf"), 8.5, 8.3),
    save_plot(legacy_power_plot(p_1to1[p_1to1$model == "make_model2", , drop = FALSE], labs, "Simulate model 2: chi-squared"),
              file.path(out$supplement, "figure-1-1to1-power-model2.pdf"), 8.5, 8.3)
  )

  d_control <- d_formula_control[d_formula_control$Method != "1:1+InMASS", , drop = FALSE]
  p_control <- p_formula_control[p_formula_control$Method != "1:1+InMASS", , drop = FALSE]
  for (form in c("ID", "MID")) {
    suffix <- if (form == "ID") "identify" else "misidentify"
    fig_no <- if (form == "ID") "2" else "3"
    d11 <- d_control[d_control$cformula == form, , drop = FALSE]
    p11 <- p_control[p_control$cformula == form, , drop = FALSE]
    method_validation <- append_plot_validation(method_validation, d11, paths, nsim, paste0("figure-", fig_no, "-controlAD-model-", suffix, "-mse"), "main")
    method_validation <- append_plot_validation(method_validation, p11, paths, nsim, paste0("figure-", fig_no, "-controlAD-model-", suffix, "-power"), "main")
    generated <- c(generated,
      save_plot(legacy_mse_plot(d11, labs), file.path(out$supplement, sprintf("figure-%s-controlAD-model-%s-mse.pdf", fig_no, suffix)), 7, 4.8),
      save_plot(legacy_power_plot(p11[p11$model == "make_model1", , drop = FALSE], labs, "Simulate model 1: normally distributed"),
                file.path(out$supplement, sprintf("figure-%s-controlAD-model-%s-power-model1.pdf", fig_no, suffix)), 8.5, 8.3),
      save_plot(legacy_power_plot(p11[p11$model == "make_model2", , drop = FALSE], labs, "Simulate model 2: chi-squared"),
                file.path(out$supplement, sprintf("figure-%s-controlAD-model-%s-power-model2.pdf", fig_no, suffix)), 8.5, 8.3)
    )
  }

  d_main <- filter_mse_power_methods(aggregate_mse_legacy(d, "Method_main"))
  p_main <- filter_mse_power_methods(aggregate_power_curve_legacy(d, "Method_main"))
  b_main <- bias_rows_legacy(d, "Method_main")
  d_main <- d_main[d_main$strata == 10, , drop = FALSE]
  p_main <- p_main[p_main$strata == 10, , drop = FALSE]
  b_k10 <- b_main[b_main$strata == 10, , drop = FALSE]
  d_main$formula <- factor(d_main$formula, levels = c(2, 1))
  p_main$formula <- factor(p_main$formula, levels = c(2, 1))
  b_k10$formula <- factor(b_k10$formula, levels = c(2, 1))
  method_validation <- append_plot_validation(method_validation, d_main, paths, nsim, "main-figure-1-mse", "main")
  method_validation <- append_plot_validation(method_validation, p_main, paths, nsim, "main-figure-2-3-power", "main")
  generated <- c(generated,
    save_plot(legacy_mse_plot(d_main, labs, stats::as.formula("formula ~ model")), file.path(out$main, "main-figure-1-mse.pdf"), 7, 4.8),
    save_plot(legacy_power_plot(p_main[p_main$model == "make_model1", , drop = FALSE], labs, "Simulate model 1: normally distributed", stats::as.formula("formula ~ n")),
              file.path(out$main, "main-figure-2-power-model1.pdf"), 9, 5.8),
    save_plot(legacy_power_plot(p_main[p_main$model == "make_model2", , drop = FALSE], labs, "Simulate model 2: chi-squared", stats::as.formula("formula ~ n")),
              file.path(out$main, "main-figure-3-power-model2.pdf"), 9, 5.8),
    save_plot(legacy_bias_plot(b_k10[b_k10$n == 40, , drop = FALSE], labs), file.path(out$main, "main-figure-0-bias.pdf"), 12, 6)
  )

  all_bias_path <- file.path(out$supplement, "main-figure-0-bias-all.pdf")
  grDevices::pdf(all_bias_path, width = 16, height = 6)
  for (st in c(5, 10, 30)) {
    for (nval in c(20, 40, 100)) {
      g4dat <- b_main[b_main$strata == st & b_main$n == nval, , drop = FALSE]
      g4dat$formula <- factor(g4dat$formula, levels = c(2, 1))
      print(legacy_bias_plot(g4dat, labs, all_bias = TRUE, title = paste(st, " studies with n=", nval, sep = "")))
    }
  }
  grDevices::dev.off()
  generated <- c(generated, all_bias_path)

  generated <- normalizePath(generated, winslash = "/", mustWork = FALSE)
  write_csv(data.frame(file = generated, stringsAsFactors = FALSE), file.path(paths$summary, sprintf("generated_figures_nsim%d.csv", nsim)))
  generated
}

make_multicov_figures <- function(paths, nsim = 10L, pilot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for figure generation.")
  }
  results_file <- file.path(paths$summary, sprintf("robustness_multicov_results_nsim%d.csv", nsim))
  if (!file.exists(results_file)) stop("Results file not found: ", results_file)

  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  d <- prepare_legacy_plot_data(results)
  d$model <- "multicov"
  d$Method <- d$Method_main
  d$formula <- factor(d$formula, levels = c(2, 1))

  labs <- legacy_labellers(d)
  labs$model <- c("Multi-covariate" = "multicov")

  out_root <- file.path(figure_paths(paths, pilot = pilot)$supplement, "robustness_multicov")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  generated <- character()

  mse_data <- filter_mse_power_methods(aggregate_mse_legacy(d, "Method_main"))
  power_data <- filter_mse_power_methods(aggregate_power_curve_legacy(d, "Method_main", b.max = 4, h = 0.1))
  bias_data <- bias_rows_legacy(d, "Method_main")
  mse_data$formula <- factor(mse_data$formula, levels = c(2, 1))
  power_data$formula <- factor(power_data$formula, levels = c(2, 1))
  bias_data$formula <- factor(bias_data$formula, levels = c(2, 1))
  power_plot_data <- power_data[power_data$strata == 10, , drop = FALSE]
  validate_mse_power_plot_methods(
    mse_data,
    file.path(paths$summary, sprintf("robustness_multicov_mse_plot_method_validation_nsim%d.csv", nsim)),
    "robustness_multicov_mse"
  )
  validate_mse_power_plot_methods(
    power_plot_data,
    file.path(paths$summary, sprintf("robustness_multicov_power_plot_method_validation_nsim%d.csv", nsim)),
    "robustness_multicov_power"
  )
  validate_power_plot_data(
    power_plot_data,
    key_cols = c("n", "formula", "Method", "sd"),
    path = file.path(paths$summary, sprintf("robustness_multicov_power_plot_diagnostics_nsim%d.csv", nsim))
  )

  generated <- c(generated,
    save_plot(
      legacy_mse_plot(mse_data, labs, stats::as.formula("formula ~ strata")),
      file.path(out_root, "robustness_multicov_mse.pdf"),
      7,
      4.8
    ),
    save_plot(
      legacy_power_plot(power_plot_data, labs, "Multi-covariate robustness scenario", stats::as.formula("formula ~ n")),
      file.path(out_root, "robustness_multicov_power.pdf"),
      9,
      5.8
    ),
    save_plot(
      legacy_bias_plot(bias_data[bias_data$n == 40 & bias_data$strata == 10, , drop = FALSE], labs),
      file.path(out_root, "robustness_multicov_bias.pdf"),
      12,
      6
    )
  )

  all_bias_path <- file.path(out_root, "robustness_multicov_bias_all.pdf")
  grDevices::pdf(all_bias_path, width = 16, height = 6)
  for (st in c(5, 10, 30)) {
    for (nval in c(20, 40, 100)) {
      g4dat <- bias_data[bias_data$strata == st & bias_data$n == nval, , drop = FALSE]
      print(legacy_bias_plot(g4dat, labs, all_bias = TRUE, title = paste(st, " studies with n=", nval, sep = "")))
    }
  }
  grDevices::dev.off()
  generated <- c(generated, all_bias_path)

  generated <- normalizePath(generated, winslash = "/", mustWork = FALSE)
  write_csv(data.frame(file = generated, stringsAsFactors = FALSE), file.path(paths$summary, sprintf("robustness_multicov_generated_figures_nsim%d.csv", nsim)))
  generated
}

make_nonlinear_figures <- function(paths, nsim = 10L, pilot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for figure generation.")
  }
  results_file <- file.path(paths$summary, sprintf("robustness_nonlinear_results_nsim%d.csv", nsim))
  if (!file.exists(results_file)) stop("Results file not found: ", results_file)

  results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  d <- prepare_legacy_plot_data(results)
  d$model <- "nonlinear"
  d$Method <- d$Method_main
  d$formula <- factor(d$formula, levels = c(2, 1))

  labs <- legacy_labellers(d)
  labs$model <- c("Nonlinear" = "nonlinear")

  out_root <- file.path(figure_paths(paths, pilot = pilot)$supplement, "robustness_nonlinear")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  generated <- character()

  mse_data <- filter_mse_power_methods(aggregate_mse_legacy(d, "Method_main"))
  power_data <- filter_mse_power_methods(aggregate_power_curve_legacy(d, "Method_main", b.max = 4, h = 0.1))
  bias_data <- bias_rows_legacy(d, "Method_main")
  mse_data$formula <- factor(mse_data$formula, levels = c(2, 1))
  power_data$formula <- factor(power_data$formula, levels = c(2, 1))
  bias_data$formula <- factor(bias_data$formula, levels = c(2, 1))
  power_plot_data <- power_data[power_data$strata == 10, , drop = FALSE]
  validate_mse_power_plot_methods(
    mse_data,
    file.path(paths$summary, sprintf("robustness_nonlinear_mse_plot_method_validation_nsim%d.csv", nsim)),
    "robustness_nonlinear_mse"
  )
  validate_mse_power_plot_methods(
    power_plot_data,
    file.path(paths$summary, sprintf("robustness_nonlinear_power_plot_method_validation_nsim%d.csv", nsim)),
    "robustness_nonlinear_power"
  )
  validate_power_plot_data(
    power_plot_data,
    key_cols = c("n", "formula", "Method", "sd"),
    path = file.path(paths$summary, sprintf("robustness_nonlinear_power_plot_diagnostics_nsim%d.csv", nsim))
  )

  generated <- c(generated,
    save_plot(
      legacy_mse_plot(mse_data, labs, stats::as.formula("formula ~ strata")),
      file.path(out_root, "robustness_nonlinear_mse.pdf"),
      7,
      4.8
    ),
    save_plot(
      legacy_power_plot(power_plot_data, labs, "Nonlinear robustness scenario", stats::as.formula("formula ~ n")),
      file.path(out_root, "robustness_nonlinear_power.pdf"),
      9,
      5.8
    ),
    save_plot(
      legacy_bias_plot(bias_data[bias_data$n == 40 & bias_data$strata == 10, , drop = FALSE], labs),
      file.path(out_root, "robustness_nonlinear_bias.pdf"),
      12,
      6
    )
  )

  all_bias_path <- file.path(out_root, "robustness_nonlinear_bias_all.pdf")
  grDevices::pdf(all_bias_path, width = 16, height = 6)
  for (st in c(5, 10, 30)) {
    for (nval in c(20, 40, 100)) {
      g4dat <- bias_data[bias_data$strata == st & bias_data$n == nval, , drop = FALSE]
      print(legacy_bias_plot(g4dat, labs, all_bias = TRUE, title = paste(st, " studies with n=", nval, sep = "")))
    }
  }
  grDevices::dev.off()
  generated <- c(generated, all_bias_path)

  generated <- normalizePath(generated, winslash = "/", mustWork = FALSE)
  write_csv(data.frame(file = generated, stringsAsFactors = FALSE), file.path(paths$summary, sprintf("robustness_nonlinear_generated_figures_nsim%d.csv", nsim)))
  generated
}

ripd_allocation_labels <- function(allocation) {
  labels <- c("1to1" = "1:1", "3to1" = "3:1", "4to0" = "Treatment only")
  unname(labels[allocation])
}

validate_ripd_truncation_plot_data <- function(plot_data, expected_files, path) {
  key_cols <- c("scenario_type", "allocation", "K", "n", "formula_id")
  duplicate_keys <- duplicated(plot_data[key_cols])
  rate_missing <- is.na(plot_data$truncation_rate)
  rate_out_of_range <- !rate_missing & (plot_data$truncation_rate < 0 | plot_data$truncation_rate > 1)
  files_created <- file.exists(expected_files)
  diagnostics <- data.frame(
    check = c(
      "truncation_rate_in_unit_interval",
      "no_duplicate_plot_rows",
      "no_missing_truncation_rate",
      "all_expected_files_created"
    ),
    passed = c(
      !any(rate_out_of_range),
      !any(duplicate_keys),
      !any(rate_missing),
      all(files_created)
    ),
    failures = c(
      sum(rate_out_of_range),
      sum(duplicate_keys),
      sum(rate_missing),
      sum(!files_created)
    ),
    stringsAsFactors = FALSE
  )
  write_csv(diagnostics, path)
  if (!all(diagnostics$passed)) {
    stop("RIPD truncation plotting validation failed. See ", path)
  }
  diagnostics
}

ripd_truncation_scenario_type <- function(simulation_family, covariate_distribution) {
  out <- rep(NA_character_, length(simulation_family))
  out[simulation_family == "main" & covariate_distribution == "normal"] <- "Normal"
  out[simulation_family == "main" & covariate_distribution == "chi2"] <- "Chi-squared"
  out[simulation_family == "robustness_multicov"] <- "Multi-covariate"
  out[simulation_family == "robustness_nonlinear"] <- "Nonlinear"
  factor(out, levels = c("Normal", "Chi-squared", "Multi-covariate", "Nonlinear"))
}

ripd_truncation_plot <- function(plot_data) {
  plot_data$allocation_label <- factor(
    ripd_allocation_labels(plot_data$allocation),
    levels = c("1:1", "3:1", "Treatment only")
  )
  plot_data$K <- factor(plot_data$K, levels = sort(unique(plot_data$K)))
  plot_data$scenario_type <- factor(
    plot_data$scenario_type,
    levels = c("Normal", "Chi-squared", "Multi-covariate", "Nonlinear")
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = truncation_rate, color = allocation_label,
                 shape = allocation_label, linetype = allocation_label, group = allocation_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    legacy_shape_scale() +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_data$n))) +
    ggplot2::scale_y_continuous(limits = c(0, NA), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Target trial sample size n",
      y = "Truncation rate",
      color = "Allocation",
      shape = "Allocation",
      linetype = "Allocation"
    ) +
    ggplot2::facet_grid(scenario_type ~ K, labeller = ggplot2::labeller(
      K = function(x) paste0(x, " studies"),
      scenario_type = ggplot2::label_value
    )) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1),
                    shape = ggplot2::guide_legend(nrow = 1),
                    fill = ggplot2::guide_legend(nrow = 1))
}

make_ripd_truncation_figures <- function(paths, nsim = 10L, pilot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for figure generation.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for truncation-rate axis labels.")
  }
  summary_file <- file.path(paths$summary, sprintf("ripd_truncation_summary_all_nsim%d.csv", nsim))
  if (!file.exists(summary_file)) stop("RIPD truncation summary file not found: ", summary_file)

  summary <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
  summary$K <- as.integer(summary$K)
  summary$n <- as.integer(summary$n)
  summary$truncation_rate <- as.numeric(summary$truncation_rate)
  summary$scenario_type <- ripd_truncation_scenario_type(summary$simulation_family, summary$covariate_distribution)
  summary <- summary[!is.na(summary$scenario_type), , drop = FALSE]
  out_root <- file.path(figure_paths(paths, pilot = pilot)$supplement, "ripd_truncation")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  specs <- data.frame(
    formula_id = c("correct", "misspecified"),
    file = c(
      "ripd_truncation_model_specified.pdf",
      "ripd_truncation_model_misspecified.pdf"
    ),
    stringsAsFactors = FALSE
  )

  old_files <- file.path(out_root, c(
    "ripd_truncation_main_model_specified.pdf",
    "ripd_truncation_main_model_misspecified.pdf",
    "ripd_truncation_multicov_model_specified.pdf",
    "ripd_truncation_multicov_model_misspecified.pdf",
    "ripd_truncation_nonlinear_model_specified.pdf",
    "ripd_truncation_nonlinear_model_misspecified.pdf"
  ))
  invisible(file.remove(old_files[file.exists(old_files)]))

  generated <- character()
  plotted <- list()
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, , drop = FALSE]
    d <- summary[
      summary$formula_id == spec$formula_id,
      ,
      drop = FALSE
    ]
    if (!nrow(d)) stop("No RIPD truncation rows for ", spec$formula_id)
    expected_types <- c("Normal", "Chi-squared", "Multi-covariate", "Nonlinear")
    missing_types <- setdiff(expected_types, as.character(unique(d$scenario_type)))
    if (length(missing_types)) {
      stop("Missing scenario types for RIPD truncation plot: ", paste(missing_types, collapse = ", "))
    }
    generated <- c(generated, save_plot(
      ripd_truncation_plot(d),
      file.path(out_root, spec$file),
      8,
      7
    ))
    plotted[[i]] <- d
  }

  generated <- normalizePath(generated, winslash = "/", mustWork = FALSE)
  plot_data <- do.call(rbind, plotted)
  validation <- validate_ripd_truncation_plot_data(
    plot_data,
    generated,
    file.path(paths$summary, sprintf("ripd_truncation_figure_validation_nsim%d.csv", nsim))
  )
  write_csv(data.frame(file = generated, stringsAsFactors = FALSE), file.path(paths$summary, sprintf("ripd_truncation_generated_figures_nsim%d.csv", nsim)))
  list(files = generated, validation = validation)
}
