## Large-data runner.
## This script expects a local CSV at data/large_fake_data.csv.

source("examples/_setup_leaveoutkss.R")

path <- file.path("data", "large_fake_data.csv")
if (!file.exists(path)) {
  stop("Place a large panel at data/large_fake_data.csv before running this example.", call. = FALSE)
}

dt <- data.table::fread(path)
data.table::setorderv(dt, cols = names(dt)[1:3])

res <- leave_out_KSS(
  y = dt[[4]],
  id = dt[[1]],
  firmid = dt[[2]],
  leave_out_level = "matches",
  type_algorithm = "JLA",
  simulations_JLA = 50,
  paral = TRUE,
  progress = TRUE
)

print_key_estimates(res)
