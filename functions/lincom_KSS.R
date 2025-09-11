lincom_KSS <- function(y, X, Z, Transform, sigma_i, labels = NULL) {

   # SET DIMENSIONS
  n <- nrow(X)
  r <- ncol(Z) + 1

  # Add Constant
  Z <- cbind(rep(1, nrow(Z)), Z)

  # COMPUTE
  xx <- t(X) %*% X # Equivalent to X'X in MATLAB
  xy <- t(X) %*% y # Equivalent to X'y in MATLAB
  beta <- sanic::solve_cg(a = xx, b = as.vector(xy), type = "CG", iter = 1000, tol = 1e-10, precond = 2, verbose = FALSE)
  wy <- Transform %*% beta
  zz <- t(Z) %*% Z # Equivalent to Z'Z in MATLAB
  any(is.na(Z))
  any(is.nan(Z))
  any(is.infinite(Z))
  any(is.na(wy))
  any(is.nan(wy))
  any(is.infinite(wy))
  numerator <- qr.solve(Z, wy)
  sigma_i <- sparseMatrix(i = 1:n, j = 1:n, x = as.vector(sigma_i),dims = c(n,n))
  error_vec <- as.vector((y - X %*% beta)^2)
  sigma_i_naive <- sparseMatrix(i = 1:n, j = 1:n, x = error_vec,dims = c(n,n))
  denominator <- numeric(r)
  denominator_naive <- numeric(r)

  for (q in 1:r) {
    v <- sparseMatrix(i = q, j = 1, x = 1, dims = c(r, 1))
    v <- qr.solve(zz, v)
    v <- Z %*% v
    v <- t(Transform) %*% v
    right <- sanic::solve_cg(a = xx, b = as.vector(v), type = "CG", iter = 1000, tol = 1e-6, precond = 2, verbose = FALSE)
    left <- t(right)
    denominator[q] <- left %*% (t(X) %*% sigma_i %*% X) %*% right
    denominator_naive[q] <- left %*% (t(X) %*% sigma_i_naive %*% X %*% right)
  }

  test_statistic <- numerator / (sqrt(denominator))

  # REPORT
  cat(strrep("*", 40), "\n")
  cat("RESULTS ON LINCOM\n")
  cat(strrep("*", 40), "\n")

  if (is.null(labels)) {
    for (q in 2:r) {
      cat(sprintf("Coefficient on - Column Number %d of Z: %f\n", q - 1, numerator[q]))
      cat(sprintf("Robust 'White' Standard Error: %f\n", sqrt(denominator_naive[q])))
      cat(sprintf("KSS Standard Error: %f\n", sqrt(denominator[q])))
      cat(sprintf("T-stat: %f\n", test_statistic[q]))
      cat(strrep("*", 40), "\n")
    }
  } else {
    for (q in 2:r) {
      tell_me <- labels[q - 1]
      cat(sprintf("Coefficient on %s: %f\n", tell_me, numerator[q]))
      cat(sprintf("Robust 'White' Standard Error: %f\n", sqrt(denominator_naive[q])))
      cat(sprintf("KSS Standard error: %f\n", sqrt(denominator[q])))
      cat(sprintf("T-stat: %f\n", test_statistic[q]))
      cat(strrep("*", 40), "\n")
    }
  }

  test_statistic <- test_statistic[-1]
  linear_combination <- numerator[-1]
  SE_linear_combination_KSS <- sqrt(denominator[-1])

  list(test_statistic = test_statistic, linear_combination = linear_combination, SE_linear_combination_KSS = SE_linear_combination_KSS)
}
