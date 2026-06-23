stable_seed <- function(..., base_seed = 1234L) {
  key <- paste(..., sep = "::")
  bytes <- utf8ToInt(key)
  value <- as.numeric(base_seed %% 2147483647L)
  for (byte in bytes) {
    value <- (value * 131 + byte) %% 2147483647
  }
  if (is.na(value) || value <= 0) value <- 1
  as.integer(value)
}

scenario_seed <- function(scenario_id, stage = "data", base_seed = 1234L) {
  stable_seed(scenario_id, stage, base_seed = base_seed)
}

replicate_seed <- function(scenario_id, replicate, stage, base_seed = 1234L) {
  stable_seed(scenario_id, replicate, stage, base_seed = base_seed)
}
