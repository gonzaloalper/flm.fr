# Test calibrated by wild bootstrap
flm.fr.test <- function(X, Y, B_boot = 500, Ad, show_plot = FALSE) {
  
  ## PREPROCESSING
  
  # Sample size
  n <- nrow(X)
  stopifnot(n == nrow(Y))
  
  #Center data - implicit column recycling
  #X <- t(t(X) - colMeans(X))
  #Y <- t(t(Y) - colMeans(Y))
  
  ## REAL WORLD
  
  # Estimate B
  B_hat <- beta_hat_net(X, Y, type = "lasso", lambda = NULL, k_folds = n/2)$beta_hat
  
  # Prediction
  Y_hat <- X %*% B_hat
  
  # Residuals
  E_hat <- Y - Y_hat
  
  # Compute statistic
  orig_stat <- PCvM_statistic(E_hat, Ad)
  
  ## BOOTSTRAP WORLD
  
  boot_stat <- numeric(B_boot)
  ones <- rep(1, n)
  
  for (i in 1:B_boot) {
    
    # Perturb residuals
    V <- fda.usc::rwild(ones, type = "golden")
    E_star <- E_hat * V # Implicit recycling, each row multiplied by the same Vi
    
    # Obtain new bootstrap observations
    Y_star <- Y_hat + E_star
    
    # Refit model - compute predictions
    B_hat_star <- beta_hat_net(X, Y_star, type = "lasso", lambda = NULL, k_folds = n/2)$beta_hat
    Y_star_hat <- X %*% B_hat_star
    
    # Residuals of refitted model
    E_star_hat <- Y_star - Y_star_hat
    
    # Compute statistic
    boot_stat[i] <- PCvM_statistic(E_star_hat, Ad)
    
  }
  
  ## P-VALUE AND PLOT
  
  # Approximation of the p-value by MC
  p_value <- mean(orig_stat < boot_stat)
  
  # Plot
  if (show_plot) {
    
    hist(boot_stat, xlim = range(c(boot_stat, orig_stat)),
         main = paste("p-value:", p_value), probability = TRUE)
    rug(boot_stat)
    abline(v = orig_stat, col = 2)
    
  }
  
  # Return
  return(list("p_value" = p_value, "orig_stat" = orig_stat, 
              "boot_stat" = boot_stat, "B_hat" = B_hat))
  
}