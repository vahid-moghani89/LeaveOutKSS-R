rsquared_comp <- function(y, id, firmid,
                        controls = NULL) {
  
  no_controls <- 0
  
  # Check for minimum required arguments
  if (missing(y) | missing(id) | missing(firmid)) {
    stop('More arguments needed')
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
  cat(strrep("*", 30), "\n", "Calculating R2\n")
  
  cat(strrep("*", 30), "\n")
  flush.console()
  tic()
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  
    DT <- data.table(y = y, id = id,old_id = id, firmid = firmid, old_firmid=firmid)
  DT_controls <- as.data.table(controls)
  
  
  
  
  # Drop stayers with a single person year observation
  DT[, count := .N, by = id]
  DT <- DT[count > 1]
  
  # Resetting ids

    DT[, firmid := match(firmid, unique(firmid))]
  
  DT[, id := match(id, unique(id))]
  DT[,interact:=paste(id,firmid,sep="-")]
  DT[,int_id:=match(interact,unique(interact))]
  
  NT <- nrow(DT)
  
  D <- sparseMatrix(i = 1:NT, j = DT$id, x = 1)
  N <- ncol(D)
  
    F <- sparseMatrix(i = 1:NT, j = DT$firmid, x = 1)
    J <- ncol(F)

  
  var_den <- var(DT$y, na.rm = TRUE)
  
  
  
  # STEP 3: Residualizing
  cat(rep("-", 25), "\n")
  cat("R2 for TWFE\n")
  cat(rep("-", 25), "\n")
  
  
    # Creating the restriction matrix S
    S <- rbind2(Diagonal(x = 1, n = J-1), Matrix::sparseMatrix(i = integer(0), j = integer(0), x=0, dims = c(1, J-1)))
    
    # Constructing the X matrix
    X <- cbind2(cbind2(D, F %*% S), as.matrix(DT_controls))
    K <- ncol(X)
  
  # Calculating xx (X' * X) and xy (X' * y)
  xx <- t(X) %*% X
  xy <- t(X) %*% DT$y
  
  
  
  # Solving for beta
  b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
  b <- as.vector(b)
  
  
  
  
  DT[, y_hat := as.vector(X %*% b)]
  
  mean_y<-mean(DT$y)
  
  DT[,residuals:=(y-y_hat)^2]
  DT[,total:=(y-mean_y)^2]
  R_squared <- 1-sum(DT$residuals)/sum(DT$total)
  
  adj_R_squared <- 1-(1-R_squared)*(NT-1)/(NT-K-1)
  cat("R2 of the TWFE Model is: ",R_squared , "\n")
  cat("Number of Parameters in the TWFE Model is: ",K , "\n")
  cat("Adjusted R2 of the TWFE Model is: ",adj_R_squared , "\n")
  
  
  
  # STEP 3: Residualizing
  cat(rep("-", 25), "\n")
  cat("R2 for Saturated model\n")
  cat(rep("-", 25), "\n")
  
  
  flag <- (DT$int_id!=1)
  D <- sparseMatrix(i = which(DT$int_id!=1), 
               j =as.integer(factor(DT[flag,]$int_id)) ,x = 1 , 
               dims=c(nrow(DT),length(unique(DT$int_id)) -1 ) )
  # Constructing the X matrix
  X <- cbind2(D, as.matrix(DT_controls))
  K <- ncol(X)
  
  # Calculating xx (X' * X) and xy (X' * y)
  xx <- t(X) %*% X
  xy <- t(X) %*% DT$y
  
  
  
  # Solving for beta
  b <- sanic::solve_cg(a = xx, b = xy, type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
  b <- as.vector(b)
  
  
  
  
  DT[, y_hat := as.vector(X %*% b)]
  
  mean_y<-mean(DT$y)
  
  DT[,residuals:=(y-y_hat)^2]
  DT[,total:=(y-mean_y)^2]
  R_squared <- 1-sum(DT$residuals)/sum(DT$total)
  
  adj_R_squared <- 1-(1-R_squared)*(NT-1)/(NT-K-1)
  cat("R2 of the Saturated Model is: ",R_squared , "\n")
  cat("Number of parameters in the Saturated Model is: ",K , "\n")
  cat("Adjusted R2 of the Saturated Model is: ",adj_R_squared , "\n")
  
}
