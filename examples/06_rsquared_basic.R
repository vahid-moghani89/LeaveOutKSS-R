## Compare TWFE and saturated-model fit on the bundled small panel.
## This example demonstrates the object returned by rsquared_comp().

source("examples/_setup_leaveoutkss.R")

dt <- read_small_test_data()
out <- example_output_stem("rsquared_basic")

res <- rsquared_comp(
  y = dt$V4,
  id = dt$V1,
  firmid = dt$V2,
  txt_file = paste0(out, ".txt"),
  progress = FALSE
)

print(res)
cat("TXT written:", file.exists(paste0(out, ".txt")), "\n")
