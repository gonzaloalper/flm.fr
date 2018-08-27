linear_model <- function(X, argvals_y, std_dev, beta)
{
  
  # Sample size
  n <- length(X)
  ly <- length(argvals_y)
  
  # Error
  noise <- matrix(rnorm(n = n * ly, mean = 0, sd = std_dev), nrow = n, ncol = ly)
  
  # Store surface
  surface_beta <- outer(X$argvals, argvals_y, beta)
  
  # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
  Y <- fdata(mdata = noise, argvals = argvals_y, rangeval = range(argvals_y), 
             names = list(main = "Functional response", ylab = "Y"))
  length_x <- diff(X$rangeval)
  
  for (j in 1:ly) {
    b <- surface_beta[, j]
    Y$data[, j] <- sapply(1:n, function(i) {
      integrateSimp1D(b * X$data[i, ], length_x)
    })
  }
  
  # Add noise
  Y$data <- Y$data + noise 
  return(Y)
  
}
