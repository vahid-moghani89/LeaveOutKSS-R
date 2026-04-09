## Compute TWFE and saturated-model R2 on the small test data.

source("examples/_setup_packages_and_functions.R")

namesrc <- file.path("data", "test.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

id     <- dt$V1
firmid <- dt$V2
year   <- dt$V3   # not used here
y      <- dt$V4

tictoc::tic()
rsquared_comp(
  y      = y,
  id     = id,
  firmid = firmid,
  filename = "rsquared_basic"
)
tictoc::toc()

## Output file:
## - rsquared_basic.txt
