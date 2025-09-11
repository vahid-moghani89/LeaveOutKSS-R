## Same as basic, but I partial out year fixed effects before KSS.

source("examples/_setup_packages_and_functions.R")

namesrc <- file.path("data", "test.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

id     <- dt$V1
firmid <- dt$V2
year   <- dt$V3
y      <- dt$V4

## year FE controls: drop the last column to avoid collinearity with intercept
controls <- model.matrix(~ factor(year) - 1)
controls <- controls[, -ncol(controls), drop = FALSE]

tictoc::tic()
res <- leave_out_KSS(
  y      = y,
  id     = id,
  firmid = firmid,
  controls = controls,
  leave_out_level = "matches",
  type_algorithm  = "JLA",
  simulations_JLA = 200,
  paral  = TRUE,
  filename = "leave_out_estimates_controls"
)
tictoc::toc()

## CSV now includes y_control_adjusted column.
