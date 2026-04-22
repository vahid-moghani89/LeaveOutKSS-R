## Small-panel lincom example using a simple region observable.
## This exercises the nested lincom result stored on the main object.

source("examples/_setup_leaveoutkss.R")

namesrc <- file.path("data", "lincom.csv")
stopifnot(file.exists(namesrc))
dt <- data.table::fread(namesrc, header = FALSE)

id     <- dt$V1
firmid <- dt$V2
# year <- dt$V3  # not used here
region <- dt$V4  # -1 for region 1, +1 for region 2
y      <- dt$V5

## region dummy: 1 if region == 1, else 0
Z_lincom <- as.numeric(region == 1)
labels_lincom <- list("Region 2 Dummy")

res <- leave_out_KSS(
  y = y,
  id = id,
  firmid = firmid,
  leave_out_level = "matches",
  type_algorithm = "JLA",
  simulations_JLA = 200,
  lincom_do = 1,
  Z_lincom = Z_lincom,
  labels_lincom = labels_lincom,
  paral = TRUE,
  progress = TRUE
)

print_key_estimates(res)
print(res$lincom)
