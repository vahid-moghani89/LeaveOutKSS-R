kss_quadratic_form <- function(sigma_i, A_1, A_2, beta, Bii) {
  # Calculate the right and left terms
  right <- A_2 %*% beta
  left <- A_1 %*% beta

  # Calculate means
  mean_left <- mean(left)
  mean_right <- mean(right)

  # Calculate the covariance manually to avoid memory issues
  theta <- sum((left - mean_left) * (right - mean_right)) / (length(left) - 1)

  # Calculate the degrees of freedom
  dof <- length(left) - 1

  # Calculate the KSS adjustment
  theta_KSS <- theta - (1 / dof) * sum(Bii * sigma_i)

  # Return the results as a list
  return(list(theta = theta, theta_KSS = theta_KSS))
}
