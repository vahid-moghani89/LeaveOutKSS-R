# leave_out_KSS: KSS leave-out variance components for two-way FE (AKM)
#
# What this does:
# - Computes KSS bias-corrected {Var(psi), Cov(alpha,psi), Var(alpha)} with optional controls.
# - Works at match-level (default) or obs-level leave-out. JLA or exact leverages.
# - Can optionally run a post-estimation lincom on firm effects vs Z.
#
# Options I care about:
#   y, id, firmid        : core inputs (vectors of same length).
#   controls             : NULL (no controls) or matrix/vec; a constant is added internally.
#   leave_out_level      : 'matches' (default, robust to within-match serial corr) or 'obs'.
#   type_algorithm       : 'JLA' (default for big n) or 'exact' (small n only).
#   simulations_JLA      : JLA sims (default 200); ↑sims = ↑accuracy, ↑time.
#   lincom_do            : 0/1. If 1, regress FE on Z (supply Z_lincom and labels).
#   Z_lincom, labels_lincom : design for lincom; labels just for nice printing.
#   filename             : CSV name for saved outputs (pe/fe + outcomes).
#   paral                : TRUE to use parallel leverages routine.
#   Cd                   : RNG seed for replicability.
#
# Notes to self:
# - Requires: data.table, Matrix, sanic (for solve_cg), and my helpers:
#   connected_set(), pruning_unbal_v3(), leverages()/leverages_parallel(),
#   kss_quadratic_form(), sigma_for_stayers(), tic()/toc().
# - I leave the verbose console prints on purpose (quick sanity checks on big runs).
leave_out_KSS <- function(y, id, firmid,
                          controls = NULL,
                          leave_out_level = 'matches',
                          type_algorithm = 'JLA',
                          simulations_JLA = 200,
                          lincom_do = 0,
                          Z_lincom = NULL,
                          labels_lincom = NULL,
                          filename = 'leave_out_estimates', paral = TRUE, Cd = 12345) {
  
  no_controls <- 0
  set.seed(Cd)
  
  # basic checks
  if (missing(y) | missing(id) | missing(firmid)) {
    stop('More arguments needed')
  }
  
  if (!leave_out_level %in% c('matches', 'obs')) {
    leave_out_level <- 'matches'
    cat("WARNING: LEAVE OUT LEVEL NOT RECOGNIZED. SWITCHED to match level.\n")
  }
  
  if (!type_algorithm %in% c('JLA', 'exact')) {
    type_algorithm <- ifelse(length(y) > 10000, 'JLA', 'exact')
    cat("WARNING: ESTIMATION TYPE NOT RECOGNIZED. SWITCHED back to default (below 10,000 observation exact and otherwise JLA).\n")
  }
  
  if (!lincom_do %in% c(0, 1)) {
    lincom_do <- 0
    cat("WARNING: lincom_do should be either 1 or 0. Switched back to default (0).\n")
  } else if (lincom_do == 1 & is.null(Z_lincom)) {
    lincom_do <- 0
    cat("WARNING: User wants to project the firm effects on some covariates Z but did not specify them; request is ignored.\n")
  }
  
  # controls handling
  if (is.vector(controls)) controls <- as.matrix(controls)
  if (is.null(controls)) {
    controls <- matrix(1, nrow=length(y), ncol=1)
    no_controls <- 1
  } else {
    controls <- cbind(1, controls)
  }
  
  # run header
  cat(strrep("*", 30), "\n", "Running KSS Correction with the following options:\n")
  switch(leave_out_level,
         matches = cat("Leave Out Strategy: Leave match out\n"),
         obs = cat("Leave Out Strategy: Leave person-year observation out\n"))
  switch(type_algorithm,
         exact = cat("Algorithm for Computation of Statistical Leverages: Exact\n"),
         JLA = cat(paste0("Algorithm for Computation of Statistical Leverages: JLA with ", simulations_JLA, " simulations.\n")))
  if (lincom_do == 1) {
    cat("Running KSS Correction after adjusting for controls as specified by the user.\n")
  }
  cat(strrep("*", 30), "\n")
  flush.console()
  tic()
  
  # data prep
  DT <- data.table(y = y, id = id, firmid = firmid)
  DT_controls <- as.data.table(controls)
  
  # STEP 1: CONNECTED SET
  DT[, lagfirmid := shift(firmid), by = id]
  DT[is.na(id) | shift(is.na(id), type = "lead", fill = TRUE), lagfirmid := NA]
  
  if (lincom_do == 1) {
    K=ncol(controls)
    DT_controls <- cbind(DT_controls, Z_lincom)
  }
  
  connected_set_result <- connected_set(y = DT$y, id = DT$id, firmid = DT$firmid, lagfirmid = DT$lagfirmid, controls = as.matrix(DT_controls))
  DT <- connected_set_result$DT
  DT_controls <- connected_set_result$DT_controls
  rm(connected_set_result); gc()
  
  # STEP 2: LEAVE-ONE-OUT CONNECTED SET
  cat(rep("-*", 14), "\n")
  cat("SECTION 1: Finding the Largest Connected Set\n")
  cat(rep("-*", 14), "\n")
  flush.console()
  
  output <- withCallingHandlers({
    pruning_unbal_v3(y = DT$y, firmid = DT$firmid, id = DT$id, id_old = DT$id_old, firmid_old = DT$firmid_old, controls = DT_controls)
  }, message = function(c) capture.output(print(c)))
  
  DT <- data.table(y = output$y,
                   firmid = output$firmid,
                   id = output$id,
                   id_old = output$id_old,
                   firmid_old = output$firmid_old)
  DT_controls <- as.data.table(output$controls)
  rm(output); gc()
  
  # drop singletons (single PY)
  DT[, count := .N, by = id]
  DT_controls <- DT_controls[DT$count > 1]
  DT <- DT[count > 1]
  
  # reindex ids/firms
  DT[, firmid := match(firmid, unique(firmid))]
  DT[, id := match(id, unique(id))]
  
  # movers/stayers flags
  DT[, first_obs := .I == 1, by = id]
  DT[, lagfirmid := shift(firmid), by = id]
  DT[first_obs==TRUE, lagfirmid := NA]
  DT[, stayer := (firmid == lagfirmid) | first_obs]
  DT[is.na(stayer), stayer:=TRUE]
  DT[, stayer_sum := sum(stayer, na.rm = TRUE), by = id]
  DT[, T := .N, by = id]
  DT[, stayer := stayer_sum == T]
  DT[, movers := !stayer]
  
  id_movers <- DT[movers==TRUE, id]
  Nmovers <- uniqueN(id_movers)
  NT <- nrow(DT)
  
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1); N <- ncol(D)
  F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1); J <- ncol(F)
  var_den <- var(DT$y, na.rm = TRUE)
  
  # summary
  cat(rep("-", 25), "\n")
  cat("Info on the leave one out connected set:\n")
  cat(rep("-", 25), "\n")
  cat("Mean of Outcome: ", mean(DT$y, na.rm = TRUE), "\n")
  cat("Variance of Outcome: ", var(DT$y, na.rm = TRUE), "\n")
  cat("# of Movers: ", Nmovers, "\n")
  cat("# of Firms: ", max(DT$firmid, na.rm = TRUE), "\n")
  cat("# of Person Year Observations: ", length(DT$y), "\n")
  cat(rep("-", 25), "\n")
  flush.console()
  
  # STEP 3: PARTIAL OUT CONTROLS
  cat(rep("-", 25), "\n")
  cat("SECTION 2: Adjusting for Controls\n")
  cat(rep("-", 25), "\n")
  
  if (lincom_do == 1) {
    Z <- as.matrix(DT_controls[, (K + 1):ncol(DT_controls), with = FALSE])
    DT_controls <- DT_controls[, 1:K, with = FALSE]
  }
  
  y_original <- DT$y
  
  if (no_controls == 0) {
    S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
    X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))
    xx <- t(X) %*% X
    xy <- t(X) %*% DT$y
    b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
    b <- as.vector(b)
    y_c_corrected <- DT$y - as.vector(X[, (N + J):ncol(X)] %*% b[(N + J):length(b)])
    DT[, y := y_c_corrected]
  }
  
  cat(rep("-", 25), "\n")
  cat("Info on the control adjusted outcomes:\n")
  cat(rep("-", 25), "\n")
  cat("Mean of Outcome: ", mean(DT$y, na.rm = TRUE), "\n")
  cat("Variance of Outcome: ", var(DT$y, na.rm = TRUE), "\n")
  
  # collapse to match means if requested
  peso <- rep(1, nrow(DT))
  y_py <- DT$y
  if (leave_out_level == 'matches') {
    DT[, match_id := as.numeric(factor(paste(id, firmid, sep = "_")))]
    DT[, peso := .N, by = match_id]
    DT <- DT[, .(
      id = mean(id),
      firmid = mean(firmid),
      id_old = mean(id_old),
      firmid_old = mean(firmid_old),
      y = mean(y),
      peso = .N
    ), by = .(match_id)]
    peso <- DT$peso
  }
  
  # STEP 4: (Pii, Bii)
  cat(rep("-", 25), "\n")
  cat("SECTION 3: COMPUTATION OF (Pii, Bii)\n")
  cat(rep("-", 25), "\n")
  flush.console()
  
  NT <- nrow(DT)
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1); N <- ncol(D)
  FM <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1); J <- ncol(FM)
  X <- cbind2(D, -FM)
  
  X_fe <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, N)) , X[, (N + 1):(N + J)])
  X_pe <- cbind2(X[, 1:N], sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, J)))
  
  row_indices_fe <- rep(1:NT, times = peso)
  row_indices_pe <- rep(1:NT, times = peso)
  X_fe <- X_fe[row_indices_fe, , drop = FALSE]
  X_pe <- X_pe[row_indices_pe, , drop = FALSE]
  
  PESO_MAT <- sparseMatrix(i = 1:NT, j = 1:NT, x = sqrt(peso))
  y_untransformed <- DT$y
  X <- PESO_MAT %*% X
  y <- PESO_MAT %*% DT$y
  xx <- t(X) %*% X
  flush.console()
  
  if (paral) {
    PBii <- leverages_parallel(X_fe = X_fe, X_pe = X_pe, X = X, xx = xx,
                               type_algorithm = type_algorithm, scale = simulations_JLA)
  } else {
    PBii <- leverages(X_fe = X_fe, X_pe = X_pe, X = X, xx = xx,
                      type_algorithm = type_algorithm, scale = simulations_JLA)
  }
  
  # STEP 5: VARIANCE COMPONENTS VIA KSS
  flush.console()
  X_fe <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, N)), FM)
  X_fe <- X_fe[rep(1:NT, times = peso), , drop = FALSE]
  X_pe <- cbind2(D, sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, J)))
  X_pe <- X_pe[rep(1:NT, times = peso), , drop = FALSE]
  S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
  X <- cbind2(D, FM %*% S)
  X <- PESO_MAT %*% X
  xx <- t(X) %*% X
  xy <- t(X) %*% y
  b <- sanic::solve_cg(a = xx, b = as.vector(xy), type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
  
  eta <- y - X %*% b
  eta_h <- eta / PBii$Mii
  sigma_i <- (y - mean(y)) * eta_h
  sigma_i <- sigma_i * PBii$correction_JLA
  
  X_fe <- X_fe[, -ncol(X_fe), drop = FALSE]
  X_pe <- X_pe[, -ncol(X_pe), drop = FALSE]
  
  if (leave_out_level == 'matches') {
    DT[, T := .N, by = id]
    DT[, stayers := T == 1]
    sigma_stayers <- sigma_stayers <- sigma_for_stayers(y = y_py,
                                                        id = DT$id,
                                                        firmid = DT$firmid,
                                                        peso = peso,
                                                        b = b)
    sigma_i[DT$stayers] <- sigma_stayers[DT$stayers]
    rm(sigma_stayers)
  }
  
  result <- kss_quadratic_form(sigma_i, X_fe, X_fe, b, PBii$Bii_fe)
  sigma_2_psi_AKM <- result[[1]]; sigma2_psi <- result[[2]]
  
  result <- kss_quadratic_form(sigma_i, X_fe, X_pe, b, PBii$Bii_cov)
  sigma_alpha_psi_AKM <- result[[1]]; sigma_psi_alpha <- result[[2]]
  
  result <- kss_quadratic_form(sigma_i, X_pe, X_pe, b, PBii$Bii_pe)
  sigma_2_alpha_AKM <- result[[1]]; sigma2_alpha <- result[[2]]
  rm(result)
  
  # STEP 6: PRINT RESULTS
  s <- "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
  cat(s, "\n"); cat(s, "\n"); cat("Decomposition Estimates:\n")
  
  s <- "PLUG-IN ESTIMATES (BIASED)"; cat(s, "\n")
  s <- paste("Variance of Firm Effects:", sigma_2_psi_AKM); cat(s, "\n")
  cat("Covariance of Firm, Person Effects:", sigma_alpha_psi_AKM, "\n")
  cat("Variance of Person Effects:", sigma_2_alpha_AKM, "\n")
  correlation <- sigma_alpha_psi_AKM / (sqrt(sigma_2_psi_AKM) * sqrt(sigma_2_alpha_AKM))
  cat("Correlation of Firm, Person Effects:", correlation, "\n")
  fraction_variance_explained <- (sigma_2_psi_AKM + 2 * sigma_alpha_psi_AKM + sigma_2_alpha_AKM) / var_den
  cat("Fraction of Variance explained by Worker and Firm Effects:", fraction_variance_explained, "\n")
  
  s <- "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
  cat(s, "\n"); cat(s, "\n")
  cat("BIAS CORRECTED ESTIMATES\n")
  cat("Variance of Firm Effects:", sigma2_psi, "\n")
  cat("Covariance of Firm and Person Effects:", sigma_psi_alpha, "\n")
  cat("Variance of Person Effects:", sigma2_alpha, "\n")
  correlation <- sigma_psi_alpha / (sqrt(sigma2_psi) * sqrt(sigma2_alpha))
  cat("Correlation of Firm, Person Effects:", correlation, "\n")
  fraction_variance_explained <- (sigma2_psi + 2 * sigma_psi_alpha + sigma2_alpha) / var_den
  cat("Fraction of Variance explained by Worker and Firm Effects:", fraction_variance_explained, "\n")
  
  # STEP 7: SAVE OUTPUT
  cat("******Creating Datasets \n")
  OrigD <- data.table(id = id, firmid = firmid)
  OrigD <- OrigD[, order_id := .I]
  setkey(OrigD,id,firmid)
  
  NewD <- data.table(id = rep(DT$id_old,times = peso),
                     pe = as.vector(X_pe %*% b),
                     firmid = rep(DT$firmid_old, times = peso),
                     fe = as.vector(X_fe %*% b))
  NewD <- unique(NewD , by = c("id","firmid"))
  setkey(NewD,id,firmid)
  
  OrigD <- NewD[OrigD, nomatch=0]
  setorder(OrigD,order_id)
  
  cat("*****Saving*****\n")
  if (no_controls==0) {
    output_dt <- data.table(
      y_original = y_original,
      y_control_adjusted = y_c_corrected,
      id = OrigD$id,
      pe = OrigD$pe,
      firm_id = OrigD$firmid,
      fe = OrigD$fe
    )
  } else {
    output_dt <- data.table(
      y_original = y_original,
      id = OrigD$id,
      pe = OrigD$pe,
      firm_id = OrigD$firmid,
      fe = OrigD$fe
    )
  }
  
  fwrite(output_dt, paste0(filename, ".csv"), row.names = FALSE, quote = FALSE)
  
  # STEP 8: LINCOM (optional)
  if (lincom_do == 1) {
    cat('Regressing the firm effects on observables...\n')
    r <- ncol(Z)
    for (q in 1:r) {
      lincom_KSS(y = y, X = X, Z = as.matrix(Z[,q]), Transform=X_fe, sigma_i = sigma_i)
    }
  }
  
  elapsed_time <- toc(log = TRUE)
}
