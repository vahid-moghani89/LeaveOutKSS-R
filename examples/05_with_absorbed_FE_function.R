## Use leave_out_KSS_fe() with a categorical control absorbed internally.
## Here I pass year as a single control column and ask the function to dummy it out.

source("examples/_setup_packages_and_functions.R")

namesrc <- file.path("data", "test.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

id     <- dt$V1
firmid <- dt$V2
year   <- dt$V3
y      <- dt$V4

controls <- cbind(year = year)

tictoc::tic()
res <- leave_out_KSS_fe(
  y      = y,
  id     = id,
  firmid = firmid,
  controls = controls,
  absorb_col = 1,
  leave_out_level = "matches",
  type_algorithm  = "JLA",
  simulations_JLA = 200,
  paral  = TRUE,
  filename = "leave_out_estimates_absorbed_fe"
)
tictoc::toc()

## Output files:
## - leave_out_estimates_absorbed_fe.csv
## - leave_out_estimates_absorbed_fe.txt
