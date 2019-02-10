linear_model <- function(X, argvals_y, std_dev, beta)
{
  
  # Sample size
  n <- dim(X$data)[1]
  ly <- length(argvals_y)
  
  # Error
  Y <- r.ou(n, t = argvals_y)
  
  # Store surface
  surface_beta <- outer(X$argvals, argvals_y, beta)
  
  # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
  #Y <- fdata(mdata = noise, argvals = argvals_y, rangeval = range(argvals_y), 
  #           names = list(main = "Functional response", ylab = "Y"))
  length_x <- diff(X$rangeval)
  
  for (j in 1:ly) {
    b <- surface_beta[, j]
    Y$data[, j] <- sapply(1:n, function(i) {
      integrateSimp1D(b * X$data[i, ], length_x)
    })
  }
  
  # Add noise
  #d <- 0
  #Y$data <- Y$data + noise + d * exp(X$data)
  return(Y)
  
}
