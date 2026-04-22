if (!requireNamespace("LeaveOutKSS", quietly = TRUE)) {
  stop(
    "Install LeaveOutKSS from CRAN before running these examples: ",
    "install.packages(\"LeaveOutKSS\")",
    call. = FALSE
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Install data.table before running these examples: ",
    "install.packages(\"data.table\")",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(LeaveOutKSS))

read_small_test_data <- function() {
  path <- system.file("extdata", "test.csv", package = "LeaveOutKSS")
  if (!nzchar(path) || !file.exists(path)) {
    stop("Could not find LeaveOutKSS example data in the installed package.", call. = FALSE)
  }

  dt <- data.table::fread(path, header = FALSE)
  data.table::setorder(dt, V1, V3)
  dt
}

print_key_estimates <- function(res) {
  print(res$estimates$table, row.names = FALSE)
  invisible(res)
}

example_output_stem <- function(stem) {
  dir_path <- file.path(tempdir(), "LeaveOutKSS_examples")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  file.path(dir_path, stem)
}
