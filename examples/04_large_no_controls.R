## Large dataset runner. No controls. Uses JLA + parallel.
## I either download a big CSV or read from ./data if already present.

source("examples/_setup_packages_and_functions.R")

## set this to an online file if you have one; otherwise, look in ./data
LARGE_URL <- "https://www.dropbox.com/s/ny5tef29ij7ran2/,!large_fake_data.csv?dl=1"
local_path <- file.path("data", "large_fake_data.csv")

if (nzchar(LARGE_URL)) {
  tmp <- tempfile(fileext = ".csv")
  message("Downloading large file..."); utils::download.file(LARGE_URL, tmp, mode = "wb")
  namesrc <- tmp
} else {
  stopifnot(file.exists(local_path))
  namesrc <- local_path
}

dt <- data.table::fread(namesrc)
## I expect columns: id, firmid, year, y (in that order)
## If names differ, adjust the indexing below.

id     <- dt[[1]]
firmid <- dt[[2]]
# year <- dt[[3]]
y      <- dt[[4]]

## JLA with fewer sims first; bump up if time allows
tictoc::tic()
res <- leave_out_KSS(
  y      = y,
  id     = id,
  firmid = firmid,
  leave_out_level = "matches",
  type_algorithm  = "JLA",
  simulations_JLA = 50,
  paral  = TRUE,
  filename = "leave_out_estimates_large"
)
tictoc::toc()

## sanity: rerun with a higher simulations_JLA and compare top-line outputs
