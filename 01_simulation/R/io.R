save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  invisible(path)
}

write_csv <- function(object, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(object, path, row.names = row.names)
  invisible(path)
}

scenario_raw_dir <- function(paths, scenario_id) {
  file.path(paths$raw, scenario_id)
}
