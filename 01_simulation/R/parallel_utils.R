make_parallel_cluster <- function(n_workers = 1L) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L) {
    stop("n_workers must be a positive integer.")
  }
  if (n_workers == 1L) return(NULL)
  cl <- parallel::makeCluster(n_workers)
  parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  cl
}

stop_parallel_cluster <- function(cluster) {
  if (!is.null(cluster)) parallel::stopCluster(cluster)
}

parallel_lapply <- function(X, FUN, ..., n_workers = 1L, cluster = NULL) {
  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1L) {
    stop("n_workers must be a positive integer.")
  }
  if (n_workers == 1L || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }
  created_cluster <- is.null(cluster)
  if (created_cluster) {
    cluster <- make_parallel_cluster(n_workers)
    on.exit(stop_parallel_cluster(cluster), add = TRUE)
  }
  parallel::parLapply(cluster, X, function(x) FUN(x, ...))
}

append_parallel_timing <- function(paths, simulation_family, nsim, n_workers, started_at, finished_at) {
  timing <- data.frame(
    simulation_family = simulation_family,
    nsim = as.integer(nsim),
    n_workers = as.integer(n_workers),
    elapsed_seconds = as.numeric(difftime(finished_at, started_at, units = "secs")),
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  path <- file.path(paths$summary, sprintf("parallel_timing_nsim%d.csv", nsim))
  if (file.exists(path)) {
    existing <- utils::read.csv(path, stringsAsFactors = FALSE)
    timing <- rbind(existing, timing)
  }
  write_csv(timing, path)
  invisible(timing)
}
