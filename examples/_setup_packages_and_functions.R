## packages I rely on; install if missing, then load
required_packages <- c(
  "Matrix", "data.table", "igraph", "tictoc", "sanic",
  "doParallel", "foreach", "slam", "spam", "dplyr"
)

to_install <- setdiff(required_packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

invisible(lapply(required_packages, require, character.only = TRUE))

shift <- data.table::shift

## source ALL function files from ./functions
fx <- list.files("functions", full.names = TRUE, pattern = "\\.R$")
stopifnot(length(fx) > 0)  # I want this to fail if the folder is empty
invisible(lapply(fx, source))

## optional: set up parallel backend for foreach (if my code uses it)
# I keep one core free
if (requireNamespace("doParallel", quietly = TRUE)) {
  nC <- max(1L, parallel::detectCores() - 1L)
  doParallel::registerDoParallel(cores = nC)
}
