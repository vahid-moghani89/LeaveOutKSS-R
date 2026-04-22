## Minimal run on the bundled small panel.
## This example demonstrates the object-returning API and opt-in progress.

source("examples/_setup_leaveoutkss.R")

dt <- read_small_test_data()

res <- leave_out_KSS(
  y = dt$V4,
  id = dt$V1,
  firmid = dt$V2,
  leave_out_level = "matches",
  type_algorithm = "JLA",
  simulations_JLA = 200,
  paral = FALSE,
  progress = TRUE
)

print(res)
print(res$sample_info$leave_one_out_connected_set)
print(utils::head(res$effects))
