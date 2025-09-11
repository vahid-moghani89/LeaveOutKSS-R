leave_out_KSS_fe <- function(y, id, firmid,
                          controls = NULL,
                          absorb_col=NULL,
                          leave_out_level = 'matches',
                          type_algorithm = 'JLA',
                          simulations_JLA = 200,
                          lincom_do = 0,
                          Z_lincom = NULL,
                          labels_lincom = NULL,
                          filename = 'leave_out_estimates', paral = TRUE,Cd = 12345) {
  
  no_controls <- 0
  set.seed(Cd)
  
  # Check for minimum required arguments
  if (missing(y) | missing(id) | missing(firmid)) {
    stop('More arguments needed')
  }
  
  # Additional checks for sizes of optional parameters
  if (!leave_out_level %in% c('matches', 'obs')) {
    leave_out_level <- 'matches'
    cat("WARNING: LEAVE OUT LEVEL NOT RECOGNIZED. SWITCHED to match level.\n")
  }
  
  if (!type_algorithm %in% c('JLA', 'exact')) {
    type_algorithm <- ifelse(length(y) > 10000, 'JLA', 'exact')
    cat("WARNING: ESTIMATION TYPE NOT RECOGNIZED. SWITCHED back to default (below 10,000 observation exact and otherwise JLA).\n")
  }
  
  # Checks related to lincom_do
  if (!lincom_do %in% c(0, 1)) {
    lincom_do <- 0
    cat("WARNING: lincom_do should be either 1 or 0. Switched back to default (0).\n")
  } else if (lincom_do == 1 & is.null(Z_lincom)) {
    lincom_do <- 0
    cat("WARNING: User wants to project the firm effects on some covariates Z but did not specify them; request is ignored.\n")
  }
  
  if (is.null(controls) & !is.null(absorb_col)){
    absorb_col = NULL
    cat("WARNING: No controls are provided, but asked for absorbing some vectors! The absorb request ignored!")
  }
  
  
  if (!is.null(controls) & !is.null(absorb_col)){
    if (!is.vector(absorb_col)) {
      absorb_col = NULL
      cat("Warning: absorv_col must be a vector, absorb command ignored!")
    }
  }
  
  if (!is.null(controls) & !is.null(absorb_col)){
    if (!all(absorb_col>=1 & absorb_col<=ncol(controls))) {
      absorb_col = NULL
      cat("Warning: absorv_col has columns out of range of columns of controls, the absorb command ignored!")
    }
  }
  
  # Check for empty controls
  if (is.vector(controls)) {
    controls <- as.matrix(controls)
  }
  
  if (is.null(controls)) {
    controls <- matrix(1, nrow=length(y), ncol=1)
    no_controls <- 1
  } else {
    controls <- cbind(1, controls)
  }
  
  # Listing options
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
  if (is.vector(absorb_col)) {
    cat("Controls include a set of categorical variables.\n")
  }
  cat(strrep("*", 30), "\n")
  flush.console()
  tic()
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  DT <- data.table(y = y, id = id, firmid = firmid)
  DT_controls <- as.data.table(controls)
  
  # STEP 1: FIND CONNECTED SET
  
  # Lagfirmid
  DT[, lagfirmid := shift(firmid), by = id]
  DT[is.na(id) | shift(is.na(id), type = "lead", fill = TRUE), lagfirmid := NA]
  
  # Update controls if lincom_do is 1
  if (lincom_do == 1) {
    K=ncol(controls)
    DT_controls <- cbind(DT_controls, Z_lincom)
  }
  
  # Find connected set
  connected_set_result <- connected_set(y = DT$y, id = DT$id, firmid = DT$firmid, lagfirmid = DT$lagfirmid, controls = as.matrix(DT_controls))
  
  # Directly updating DT and DT_controls
  DT <- connected_set_result$DT
  DT_controls <- connected_set_result$DT_controls
  rm(connected_set_result)
  gc()
  
  # STEP 2: LEAVE ONE OUT CONNECTED SET
  # Display section header
  cat(rep("-*", 14), "\n")
  cat("SECTION 1: Finding the Largest Connected Set\n")
  cat(rep("-*", 14), "\n")
  flush.console()
  # Invoke the pruning_unbal_v3 function and update datasets
  # Invoke the pruning_unbal_v3 function
  output <- withCallingHandlers({
    pruning_unbal_v3(y = DT$y, firmid = DT$firmid, id = DT$id, id_old = DT$id_old, firmid_old = DT$firmid_old, controls = DT_controls)
  }, message = function(c) capture.output(print(c)))
  
  # Create new data.tables
  DT <- data.table(y = output$y,
                   firmid = output$firmid,
                   id = output$id,
                   id_old = output$id_old,
                   firmid_old = output$firmid_old)
  
  DT_controls <- as.data.table(output$controls)
  rm(output)
  gc()
  
  # Drop stayers with a single person year observation
  DT[, count := .N, by = id]
  DT_controls<-DT_controls[DT$count > 1]
  DT <- DT[count > 1]
  
  
  # Resetting ids
  DT[, firmid := match(firmid, unique(firmid))]
  DT[, id := match(id, unique(id))]
  
  # Important Auxiliaries
  # Create a flag for the first observation for each worker
  DT[, first_obs := .I == 1, by = id]
  
  # Create a lagged firmid, with NA for the first observation per worker
  DT[, lagfirmid := shift(firmid), by = id]
  DT[first_obs==TRUE, lagfirmid := NA]
  
  # Identify stayers: same firmid as lagfirmid, flag first observation as stayer
  DT[, stayer := (firmid == lagfirmid) | first_obs]
  DT[is.na(stayer), stayer:=TRUE]
  
  # Sum up stayers by id to determine if all records are stayers
  DT[, stayer_sum := sum(stayer, na.rm = TRUE), by = id]
  DT[, T := .N, by = id]
  DT[, stayer := stayer_sum == T]
  
  # Identify movers: not a stayer in any record
  DT[, movers := !stayer]
  
  # Calculate the number of unique movers
  id_movers <- DT[movers==TRUE, id]
  Nmovers <- uniqueN(id_movers) # Number of unique movers
  
  NT <- nrow(DT)
  
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
  J <- ncol(F)
  var_den <- var(DT$y, na.rm = TRUE)
  
  # Summarize
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
  
  
  # STEP 3: Residualizing and Collapsing
  cat(rep("-", 25), "\n")
  cat("SECTION 2: Adjusting for Controls\n")
  cat(rep("-", 25), "\n")
  
  #cat("Number of rows of DT:", nrow(DT),"\n")
  #cat("Number of rows of DT_controls:", nrow(DT_controls),"\n")
  # Assuming D, F, controls, y, J, N, and no_controls are already defined
  if (lincom_do == 1) {
    Z <- as.matrix(DT_controls[, (K + 1):ncol(DT_controls), with = FALSE])
    DT_controls <- DT_controls[, 1:K, with = FALSE]
  }
  
  y_original <- DT$y
  
  if (no_controls == 0) {
    # Creating the restriction matrix S
    S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
    if (is.vector(absorb_col)){
      absorb_col = absorb_col +1 #to adjust for the fact that a vector of 1 is added as the first column
      }
    
    if (is.vector(absorb_col)){
      cat_controls<-DT_controls[,..absorb_col,with=FALSE]
      DT_controls<-as.matrix(DT_controls[,-absorb_col,with=FALSE])
      
      for (i in seq_len(ncol(cat_controls))) {
        cat_controls[,i]<-match(cat_controls[[i]],unique(cat_controls[[i]]))
        flag <- (cat_controls[[i]]!=1)
        dummy_control <- sparseMatrix(i = which(cat_controls[[i]]!=1), 
                                      j =as.integer(factor(cat_controls[flag,][[i]])) ,x = 1 , 
                                      dims=c(nrow(cat_controls),length(unique(cat_controls[[i]])) -1 ) )
        DT_controls <-cbind2(DT_controls,dummy_control)
        rm(dummy_control)
        gc()
      } 
      cat("Categorical variables turned into dummies successfully!\n")
      flush.console()
    } else {
      DT_controls<-as.matrix(DT_controls)
    }
    # Constructing the X matrix
    X <- cbind2(cbind2(D, F %*% S), DT_controls)
    
    # Calculating xx (X' * X) and xy (X' * y)
    xx <- t(X) %*% X
    xy <- t(X) %*% DT$y
    
    
    
    # Solving for beta
    b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = TRUE)
    b <- as.vector(b)
    
    
    
    # Variance decomposition based on residualized outcome
    y_c_corrected<-DT$y - as.vector(X[, (N + J):ncol(X)] %*% b[(N + J):length(b)])
    DT[, y := y_c_corrected]
  }
  
  
  # Summarize after adjusting for controls
  cat(rep("-", 25), "\n")
  cat("Info on the control adjusted outcomes:\n")
  cat(rep("-", 25), "\n")
  cat("Mean of Outcome: ", mean(DT$y, na.rm = TRUE), "\n")
  cat("Variance of Outcome: ", var(DT$y, na.rm = TRUE), "\n")
  
  # Collapsing
  peso <- rep(1, nrow(DT))
  y_py <- DT$y
  
  if (leave_out_level == 'matches') {
    # Generate match_id
    DT[, match_id := as.numeric(factor(paste(id, firmid, sep = "_")))]
    
    # Calculate weights based on the length of a given match
    DT[, peso := .N, by = match_id]
    
    # Collapse data down to match-means
    DT <- DT[, .(
      id = mean(id),
      firmid = mean(firmid),
      id_old = mean(id_old),
      firmid_old = mean(firmid_old),
      y = mean(y),
      peso = .N  # This counts the number of rows in each group
    ), by = .(match_id)]
    peso <- DT$peso
  }
  
  # STEP 4: COMPUTATION OF (Pii, Bii)
  cat(rep("-", 25), "\n")
  cat("SECTION 3: COMPUTATION OF (Pii, Bii)\n")
  cat(rep("-", 25), "\n")
  flush.console()
  
  # Build Design
  NT <- nrow(DT)
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  FM <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
  J <- ncol(FM)
  X <- cbind2(D, -FM)
  
  # Debugging messages
  # Weighting Matrices
  X_fe <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, N)) , X[, (N + 1):(N + J)])
  X_pe <- cbind2(X[, 1:N], sparseMatrix(i = integer(0), j = integer(0), x=0 , dims = c(NT, J)))
  # Replicate rows according to the weights 'peso'
  row_indices_fe <- rep(1:NT, times = peso)
  row_indices_pe <- rep(1:NT, times = peso)
  
  X_fe <- X_fe[row_indices_fe, , drop = FALSE]
  X_pe <- X_pe[row_indices_pe, , drop = FALSE]
  
  # FGLS Transformation
  PESO_MAT <- sparseMatrix(i = 1:NT, j = 1:NT, x = sqrt(peso))
  
  
  y_untransformed <- DT$y
  
  X <- PESO_MAT %*% X
  
  y <- PESO_MAT %*% DT$y
  
  # Compute xx (X' * X)
  xx <- t(X) %*% X
  flush.console()
  if (paral) {
    PBii <- leverages_parallel(X_fe = X_fe, X_pe = X_pe, X = X, xx = xx,
                               type_algorithm = type_algorithm, scale = simulations_JLA)
  } else {
    PBii <- leverages(X_fe = X_fe, X_pe = X_pe, X = X, xx = xx,
                      type_algorithm = type_algorithm, scale = simulations_JLA)
  }
  
  
  
  # Extract the results
  #Pii <- results$Pii
  #Mii <- results$Mii
  #correction_JLA <- results$correction_JLA
  #Bii_fe <- results$Bii_fe
  #Bii_cov <- results$Bii_cov
  #Bii_pe <- results$Bii_pe
  
  # STEP 5: ESTIMATION OF VARIANCE COMPONENTS
  # We use the statistical leverages, Pii, and the Bii associated with a given
  # variance component to bias correct these quantities using the KSS approach
  flush.console()
  # Bind feature matrices and add debugging message
  X_fe <- cbind2(sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, N)), FM)
  
  X_fe <- X_fe[rep(1:NT, times = peso), , drop = FALSE]
  
  
  # Bind prediction matrices and add debugging message
  X_pe <- cbind2(D, sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(NT, J)))
  X_pe <- X_pe[rep(1:NT, times = peso), , drop = FALSE]
  
  # Bind S matrix and add debugging message
  S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
  
  # Combine D and F matrices, and add debugging message
  X <- cbind2(D, FM %*% S)
  X <- PESO_MAT %*% X
  
  # Compute crossproducts and add debugging messages
  xx <- t(X) %*% X
  xy <- t(X) %*% y
  
  # Execute the Conjugate Gradient solver and add debugging message
  b <- sanic::solve_cg(a = xx, b = as.vector(xy), type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)

  
  # Residual calculation and KSS estimate
  eta <- y - X %*% b
  eta_h <- eta / PBii$Mii  # Leave one out residual
  sigma_i <- (y - mean(y)) * eta_h  # KSS estimate of individual variance
  sigma_i <- sigma_i * PBii$correction_JLA # Adjustment for non-linear bias
  
  # Drop the last column from the fixed and random effects design matrices
  X_fe <- X_fe[, -ncol(X_fe), drop = FALSE]
  X_pe <- X_pe[, -ncol(X_pe), drop = FALSE]
  
  if (leave_out_level == 'matches') {
    # Calculate the number of times each 'id' is present in the dataset
    DT[, T := .N, by = id]
    
    
    # Determine 'stayers' (where T == 1)
    DT[, stayers := T == 1]
    # Calculate 'sigma_stayers' using your 'sigma_for_stayers' function
    # Assuming 'sigma_for_stayers' is a vectorized function that returns a vector of the same length
    sigma_stayers <- sigma_stayers <- sigma_for_stayers(y = y_py,
                                                        id = DT$id,
                                                        firmid = DT$firmid,
                                                        peso = peso,
                                                        b = b)
    sigma_i[DT$stayers]<-sigma_stayers[DT$stayers]
    rm(sigma_stayers)
    # Update 'sigma_i' for stayers
  }
  
  
  
  

  
  # Compute the quadratic form for the first set and assign the values to variables
  result <- kss_quadratic_form(sigma_i, X_fe, X_fe, b, PBii$Bii_fe)
  sigma_2_psi_AKM <- result[[1]]
  sigma2_psi <- result[[2]]
  # Compute the quadratic form for the second set and assign the values to variables
  result <- kss_quadratic_form(sigma_i, X_fe, X_pe, b, PBii$Bii_cov)
  sigma_alpha_psi_AKM <- result[[1]]
  sigma_psi_alpha <- result[[2]]
  
  # Compute the quadratic form for the third set and assign the values to variables
  result <- kss_quadratic_form(sigma_i, X_pe, X_pe, b, PBii$Bii_pe)
  sigma_2_alpha_AKM <- result[[1]]
  sigma2_alpha <- result[[2]]
  rm(result) # Remove the result from memory
  
  # STEP 6: PRINTING RESULTS
  s <- "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
  cat(s, "\n")
  cat(s, "\n")
  cat("Decomposition Estimates:\n")
  
  s <- "PLUG-IN ESTIMATES (BIASED)"
  cat(s, "\n")
  s <- paste("Variance of Firm Effects:", sigma_2_psi_AKM)
  cat(s, "\n")
  cat("Covariance of Firm, Person Effects:", sigma_alpha_psi_AKM, "\n")
  
  cat("Variance of Person Effects:", sigma_2_alpha_AKM, "\n")
  correlation <- sigma_alpha_psi_AKM / (sqrt(sigma_2_psi_AKM) * sqrt(sigma_2_alpha_AKM))
  cat("Correlation of Firm, Person Effects:", correlation, "\n")
  fraction_variance_explained <- (sigma_2_psi_AKM + 2 * sigma_alpha_psi_AKM + sigma_2_alpha_AKM) / var_den
  cat("Fraction of Variance explained by Worker and Firm Effects:", fraction_variance_explained, "\n")
  
  s <- "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
  cat(s, "\n")
  cat(s, "\n")
  cat("BIAS CORRECTED ESTIMATES\n")
  cat("Variance of Firm Effects:", sigma2_psi, "\n")
  cat("Covariance of Firm and Person Effects:", sigma_psi_alpha, "\n")
  
  cat("Variance of Person Effects:", sigma2_alpha, "\n")
  correlation <- sigma_psi_alpha / (sqrt(sigma2_psi) * sqrt(sigma2_alpha))
  cat("Correlation of Firm, Person Effects:", correlation, "\n")
  fraction_variance_explained <- (sigma2_psi + 2 * sigma_psi_alpha + sigma2_alpha) / var_den
  cat("Fraction of Variance explained by Worker and Firm Effects:", fraction_variance_explained, "\n")
  
  

  # STEP 7: SAVING OUTPUT
  # Combine data into a data.table
  cat("******Creating Datasets \n")
  OrigD <- data.table(id = id, firmid = firmid)
  OrigD <- OrigD[, order_id := .I]
  setkey(OrigD,id,firmid)
  
  NewD <- data.table(id = rep(DT$id_old,times = peso),
                     pe = as.vector(X_pe %*% b),  # Ensuring the result is a vector
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
      pe = OrigD$pe,  # Ensuring the result is a vector
      firm_id = OrigD$firmid,
      fe = OrigD$fe
    )
  } else{
    output_dt <- data.table(
      y_original = y_original,
      id = OrigD$id,
      pe = OrigD$pe,  # Ensuring the result is a vector
      firm_id = OrigD$firmid,
      fe = OrigD$fe
    )
  }
  
  
  # Save the data.table as a CSV file
  fwrite(output_dt, paste0(filename, ".csv"), row.names = FALSE, quote = FALSE)
  
  
  
  # STEP 8: LINCOM
  if (lincom_do == 1) {
    cat('Regressing the firm effects on observables...\n')
    
    # Check if labels_lincom is not NULL
    r <- ncol(Z)
    for (q in 1:r) {
      lincom_KSS(y = y, X = X, Z = as.matrix(Z[,q]), Transform=X_fe, sigma_i = sigma_i)  
    }
    
  }
  
  
  
  elapsed_time <- toc(log = TRUE)
}
