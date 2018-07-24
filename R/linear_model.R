linear_model = function(X, n, ly, argvals_y, std_dev, surface_beta)
{
  # Error
  noise <- matrix(rnorm(n = n * ly, mean = 0, sd = std_dev), nrow = n, ncol = ly)
  
  # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
  obj <- fdata(mdata = noise, argvals = argvals_y)
  
  #Y$data <- X$data %*% surface_beta * hx
  
  for (j in 1:ly){
    b = surface_beta[,j]
    sapply(1:n, function(i){
    obj[i,j] <<- sdetorus::integrateSimp1D(b * X$data[i,], X_range)
    })
  }
  
  Y$data <- Y$data + noise 
  rangevals_y <- c(tt[1], tail(tt, n=1))
  
  Y = fdata(data = obj$data, argvals = tt, rangeval = rangevals_y, 
            names = list(main = "Functional response", ylab = Y))
  return(Y)
}
