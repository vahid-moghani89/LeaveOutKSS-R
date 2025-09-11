## After estimating FE, regress FE on a Z (e.g., region dummy).
## This follows the layout of lincom.csv (V1..V5).

source("examples/_setup_packages_and_functions.R")

namesrc <- file.path("data", "lincom.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

id     <- dt$V1
firmid <- dt$V2
# year <- dt$V3  # not using here
region <- dt$V4  # -1 for region 1, +1 for region 2
y      <- dt$V5

## region dummy: 1 if region==1, else 0
Z_lincom <- as.numeric(region == 1)
labels_lincom <- list("Region 2 Dummy")

tictoc::tic()
res <- leave_out_KSS(
  y      = y,
  id     = id,
  firmid = firmid,
  leave_out_level = "matches",
  type_algorithm  = "JLA",
  simulations_JLA = 100,
  lincom_do = 1,
  Z_lincom  = Z_lincom,
  labels_lincom = labels_lincom,
  paral  = TRUE,
  filename = "leave_out_estimates_lincom"
)
tictoc::toc()

## console should show lincom output after the main decomposition
