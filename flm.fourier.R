flm.fourier <- function(B)
{
  ## FLM WITH FUNCTIONAL RESPONSE EXAMPLE
  
  library(viridis)
  library(fda.usc)
  
  # Execute all the functions in the R directory from the above file:
  #sourceDir <- function(directory = "R", ...) {
  #  invisible(sapply(dir(directory),
  #                   function(x) source(paste0(directory, "/", x), ...)))
  #}
  #sourceDir()
  
  # Grids of X and Y
  n = 100
  lx = 201; argvalsx <- seq(0, 1, l = lx) # Argvals of X
  ly = 201; argvalsy <- seq(0, 1, l = ly) # Argvals of Y
  
  # Data generation
  data <- data_generation(n = n, ss = argvalsx, tt = argvalsy)
  X <- data[[1]]; Y <- data[[2]]; beta <- data[[3]]; surface_beta <- data[[4]]
  
  # Plots
  par(mfrow = c(1, 2))
  plot(X)
  plot(Y)
  
  # Basis expansions
  kx <- 15; ky <- 15
  fourier <- fourier_expansion(X, Y, kx, ky)
  projX <- fourier[[1]]; projY <- fourier[[2]]
  basis_fourier_x <- fourier[[3]]; basis_fourier_y <- fourier[[4]]
  
  # Coefficients (estimated and theoretical)
  hat_b <- pseudoinverse(t(projX) %*% projX) %*% t(projX) %*% projY
  coefs_beta <- matrix(0, nrow = kx, ncol = ky)
  for (i in 1:kx){
    for (j in 1:ky){
      outer_product <- outer(basis_fourier_x[,i],basis_fourier_y[,j])
      acc <- surface_beta * outer_product
      coefs_beta[i,j] <- integrateSimp2D(acc, c(X_range, Y_range))
    }
  }
  
  integrateSimp2D(outer(basis_fourier_x[,i],basis_fourier_y[,j])^2,
                            c(X_range,Y_range))
  
  surface_fourier_hat <- matrix(0, nrow = lx, ncol = ly)
  surface_fourier <- matrix(0, nrow = lx, ncol = ly)
  
  for (i in 1:kx){
    for (j in 1:ky){
      outer_product <- outer(basis_fourier_x[,i],basis_fourier_y[,j])
      acc1 <- coefs_beta[i,j] * outer_product
      surface_fourier <- surface_fourier + acc1
      acc2 <- hat_b[i,j] * outer_product
      surface_fourier_hat <- surface_fourier_hat + acc2
    }
  }
  
  par(mfrow = c(1, 3))
  minima <- c(min(surface_fourier), min(surface_beta), min(surface_fourier_hat))
  maxima <- c(max(surface_fourier), max(surface_beta), max(surface_fourier_hat))
  limits = c(min(minima), max(maxima))
  image(ss, tt, surface_beta, col = viridis(20), zlim = limits)
  image(ss, tt, surface_fourier, col = viridis(20), zlim = limits)
  image(ss, tt, surface_fourier_hat, col = viridis(20), zlim = limits)
  
  ## RESIDUALS AND BOOTSTRAP
  y_hat <- Y; y_hat2 <- Y
  #y_hat$data <- projX %*% hat_b
  y_hat2$data <- projY%*%t(basis_fourier_y)
  y_hat$data <- projX %*% hat_b %*% t(basis_fourier_y)
  par(mfrow = c(1, 4))
  plot(Y)
  plot(y_hat2)
  plot(y_hat)
  fresiduals <- y_hat - Y
  plot(fresiduals)
  
  res_star <- fresiduals; for (j in 1:dim(fresiduals$data)[1]){ #PERTURBED RESIDUALS
    res_star$data[j,] <- fresiduals$data[j,] * rwild(1, "golden")}
  Y_star <- Y - fresiduals + res_star
  
  ## WILD BOOTSTRAP   RESAMPLING
  
  P = pseudoinverse(t(projX) %*% projX) %*% t(projX)
  PCvM_star <- 0; 
  
  # Adot
  Adot.vec=Adot(X)
  
  # Obtain the entire matrix Adot
  Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
  Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
  Ad=t(Ad)
  Ad=Ad+t(Ad)-diag(diag(Ad))
  
  PCvM <- PCvM_statistic(projX, fresiduals$data %*% basis_fourier_y, Ad)
  
  res_star <- fresiduals
  
  Y_star <- Y; y_hat_star <- Y; res_hat_star <- fresiduals
  projY_star <- projY
  
  for (i in 1:B) {
    for (j in 1:dim(fresiduals$data)[1]){
      res_star$data[j,] <- fresiduals$data[j,] * rwild(1, "golden")
    }
    
    Y_star <- Y - fresiduals + res_star
    
    for (i in 1:n){
      for (j in 1:ky){
        projY_star[i,j] <- integrateSimp1D(Y_star$data[i,]*basis_fourier_y[,j],Y_range)
      }
    }
    
    hat_b_star <- P %*% projY_star
    
    y_hat_star$data <- projX %*% hat_b_star %*% t(basis_fourier_y)
    
    res_hat_star <- Y_star - y_hat_star
    PCvM_star[i] = na.omit(PCvM_statistic(pcX$x, res_hat_star$data %*% basis_fourier_y, Ad))
  }
  
  PCvM_star <- na.omit(PCvM_star)
  par(mfrow = c(1,1))
  hist(PCvM_star)
  abline(v = PCvM, col = 2)
  pvalue = sum(PCvM_star > PCvM)/length(PCvM_star); print(pvalue)
  
  par(mfrow = c(2, 4))
  plot(Y)
  plot(y_hat)
  plot(Y_star)
  plot(y_hat_star)
  plot(fresiduals)
  plot(fresiduals)
  plot(res_star)
  plot(res_hat_star)
  
  return(pvalue)
}
