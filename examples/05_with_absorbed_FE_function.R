## Small-panel run using leave_out_KSS_fe() to absorb year effects internally.
## This exercises the absorbed-control pathway without pre-building dummies.

source("examples/_setup_leaveoutkss.R")

dt <- read_small_test_data()

res <- leave_out_KSS_fe(
  y = dt$V4,
  id = dt$V1,
  firmid = dt$V2,
  controls = cbind(year = dt$V3),
  absorb_col = 1,
  leave_out_level = "matches",
  type_algorithm = "JLA",
  simulations_JLA = 200,
  paral = FALSE,
  progress = TRUE
)

print_key_estimates(res)
print(utils::head(res$effects))
