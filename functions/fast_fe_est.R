fast_fe_est <- function(y, id, firmid = NULL,
                          controls = NULL,
                          filename = 'fe_est.csv') {
  
  no_controls <- 0
  
  # Check for minimum required arguments
  if (missing(y) | missing(id)) {
    stop('More arguments needed')
  }
  if (is.null(firmid)){
    cat('running a one way fixed effect model.\n')
  } else {
    cat('running a two way fixed effect model.\n')
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
  cat(strrep("*", 30), "\n", "Running a Fixed Effect Model\n")
  
  cat(strrep("*", 30), "\n")
  flush.console()
  tic()
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  if (is.null(firmid)) {
    DT <- data.table(y = y, id = id, old_id=id)
  } else {
    DT <- data.table(y = y, id = id,old_id = id, firmid = firmid, old_firmid=firmid)
  }
  DT_controls <- as.data.table(controls)
  

  
  
  # Drop stayers with a single person year observation
  DT[, count := .N, by = id]
  DT <- DT[count > 1]
  
  # Resetting ids
  if (!is.null(firmid)){
    DT[, firmid := match(firmid, unique(firmid))]
  }
  DT[, id := match(id, unique(id))]
  
  NT <- nrow(DT)
  
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  if (!is.null(firmid)){
    F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
    J <- ncol(F)
  }
 
  var_den <- var(DT$y, na.rm = TRUE)
  
  
  
  # STEP 3: Residualizing
  cat(rep("-", 25), "\n")
  cat("Residualizing\n")
  cat(rep("-", 25), "\n")
  
  
  if (!is.null(firmid)){
    # Creating the restriction matrix S
    S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
    
    # Constructing the X matrix
    X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))
  } else {
    X <- cbind2(D, as.matrix(DT_controls))
  }
    # Calculating xx (X' * X) and xy (X' * y)
    xx <- t(X) %*% X
    xy <- t(X) %*% DT$y
    
    
    
    # Solving for beta
    b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
    b <- as.vector(b)
    
    
    
  
    DT[, y_adj := y - as.vector(X[, (N + J):ncol(X)] %*% b[(N + J):length(b)])]
    DT[, y_hat := as.vector(X %*% b)]
  # SAVING OUTPUT
  # Combine data into a data.table
  cat("Saving\n")
  if (!is.null(firmid)){
  output_dt <- data.table(
    y_hat = DT$y_hat,
    y_adj = DT$y_adj,
    personid = DT$old_id,
    firm_id = DT$old_firmid)
  } else {
    output_dt <- data.table(
      y_hat = DT$y_hat,
      y_adj = DT$y_adj,
      personid = DT$old_id)
  }
  
  # Save the data.table as a CSV file
  fwrite(output_dt, paste0(filename), row.names = FALSE, quote = FALSE)
  
}
