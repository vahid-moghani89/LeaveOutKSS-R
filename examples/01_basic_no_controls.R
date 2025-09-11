## Minimal example: no controls, match-level leave-out, JLA
## I source functions first, then run. Keep it simple.

source("examples/_setup_packages_and_functions.R")

## read small test data
namesrc <- file.path("data", "test.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

## columns (stick to the MATLAB test layout)
id     <- dt$V1
firmid <- dt$V2
year   <- dt$V3   # not used here
y      <- dt$V4

## run KSS
tictoc::tic()
res <- leave_out_KSS(
  y      = y,
  id     = id,
  firmid = firmid,
  leave_out_level = "matches",
  type_algorithm  = "JLA",
  simulations_JLA = 200,
  paral  = TRUE,
  filename = "leave_out_estimates_basic"
)
tictoc::toc()

## notes to self:
## - printed output has both plug-in and bias-corrected pieces
## - a CSV "leave_out_estimates_basic.csv" is created in working dir
