## Small-panel run with year fixed effects partialled out through a control matrix.
## This example also demonstrates opt-in file export to tempdir().

source("examples/_setup_leaveoutkss.R")

dt <- read_small_test_data()

controls <- model.matrix(~ factor(dt$V3) - 1)
controls <- controls[, -ncol(controls), drop = FALSE]
out <- example_output_stem("leave_out_controls")

res <- leave_out_KSS(
  y = dt$V4,
  id = dt$V1,
  firmid = dt$V2,
  controls = controls,
  leave_out_level = "matches",
  type_algorithm = "JLA",
  simulations_JLA = 200,
  paral = FALSE,
  csv_file = paste0(out, ".csv"),
  txt_file = paste0(out, ".txt"),
  progress = FALSE
)

print_key_estimates(res)
cat("CSV written:", file.exists(paste0(out, ".csv")), "\n")
cat("TXT written:", file.exists(paste0(out, ".txt")), "\n")
