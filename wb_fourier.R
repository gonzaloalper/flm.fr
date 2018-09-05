library(viridis)
library(fda.usc)

#Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Proof-of-concept illustration of Wild bootstrap resampling in multivariate 
# linear regression with multivariate response

# Simulate data from model

# Grids of X and Y
n = 100
lx = 201; argvalsx <- seq(0, 2, l = lx) # Argvals of X
ly = 201; argvalsy <- seq(0, 3, l = ly) # Argvals of Y

# Data generation
X <- rproc2fdata(n = n, argvalsx, sigma="brownian")

# Theoretical beta
beta <- function(s, t) {
  cos(2 * pi * t * s) / (1 + (s - 0.5)^2)}

# Visualization
surface_beta <- outer(argvalsx, argvalsy, FUN = beta)

# Error
noise <- matrix(rnorm(n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = noise, argvals = argvalsy, rangeval = range(argvalsy), 
           names = list(main = "Functional response", ylab = "Y"))
length_x <- diff(X$rangeval)

for (j in 1:ly) {
  b <- surface_beta[, j]
  Y$data[, j] <- sapply(1:n, function(i) {
    integrateSimp1D(b * X$data[i, ], length_x)
  })
}

# Add noise
d <- 0.05
Y1 <- Y
Y1$data <- Y$data + noise + d * exp(X$data)

# # Plots
# X$data <- X$data[1:10,]
# Y$data <- Y$data[1:10,]
par(mfrow = c(1, 2))
plot(X,xlab = "s", ylab = "X(s)",ylim=c(-5,5))
plot(Y1,xlab = "t", ylab = "Y(t)",ylim=c(-2,2))

# Basis expansions
kx <- 7; ky <- 7
X_range <- X$rangeval[2] - X$rangeval[1]
Y_range <- Y$rangeval[2] - Y$rangeval[1]

basis_fourier_x <- t(eval.basis(evalarg = argvalsx,
                                basisobj = create.fourier.basis(rangeval = c(X$rangeval[1], X$rangeval[2]), 
                                                                nbasis = kx)))
basis_fourier_y <- t(eval.basis(evalarg = argvalsy,
                                basisobj = create.fourier.basis(rangeval = c(Y$rangeval[1], Y$rangeval[2]), 
                                                                nbasis = ky)))

for (i in 1:kx){
  basis_fourier_x[i,] <-
    basis_fourier_x[i,]/sqrt(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
}

for (i in 1:ky){
  basis_fourier_y[i,] <-
    basis_fourier_y[i,]/sqrt(integrateSimp1D(basis_fourier_y[i,]^2, c(Y_range)))
}

basis_fourier_x <- t(basis_fourier_x); basis_fourier_y <- t(basis_fourier_y)

projX <- matrix(0, nrow = n, ncol = kx)
projY <- matrix(0, nrow = n, ncol = ky)

for (i in 1:n){
  for (j in 1:kx){
    projX[i,j] <- integrateSimp1D(X$data[i,]*basis_fourier_x[,j],X_range)
  }
}

for (i in 1:n){
  for (j in 1:ky){
    projY[i,j] <- integrateSimp1D(Y$data[i,]*basis_fourier_y[,j],Y_range)
  }
}

# Adot
Adot.vec=Adot(X)

# Obtain the entire matrix Adot
Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
Ad=t(Ad)
Ad=Ad+t(Ad)-diag(diag(Ad))

# Test calibrated by wild bootstrap
test <- function(X, Y, B_boot = 500, Ad, show_plot = FALSE) {
  
  ## PREPROCESSING
  
  # Sample size
  n <- nrow(X)
  stopifnot(n == nrow(Y))
  
  #Center data - implicit column recycling
  #X <- t(t(X) - colMeans(X))
  #Y <- t(t(Y) - colMeans(Y))
  
  ## REAL WORLD
  
  # Estimate B
  H <- pseudoinverse(crossprod(X)) %*% t(X) 
  B_hat <- H %*% Y
  
  # Hat matrix
  H <- X %*% H
  
  # Prediction
  Y_hat <- H %*% Y

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
    
    # Refit model - compute predictions directly
    Y_star_hat <- H %*% Y_star 
    
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

# Individual check
test(projX, projY, B_boot <-  500, Ad)$p_value

# Monte Carlo
M <- 500
t <- list()
pb <- txtProgressBar(style = 3)
for (counter in 1:M) {
  # Data generation
  n <- 100
  kx <- 7; ky <- 7
  X <- rproc2fdata(n = n, argvalsx, sigma="brownian")
  
  # Error
  noise <- matrix(rnorm(n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)
  
  # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
  Y <- fdata(mdata = noise, argvals = argvalsy, rangeval = range(argvalsy), 
             names = list(main = "Functional response", ylab = "Y"))
  length_x <- diff(X$rangeval)
  
  for (j in 1:ly) {
    b <- surface_beta[, j]
    Y$data[, j] <- sapply(1:n, function(i) {
      integrateSimp1D(b * X$data[i, ], length_x)
    })
  }
  
  # Add noise
  d <- 0.002
  Y$data <- Y$data + noise + d * exp(X$data)
  
  # # Plots
  # par(mfrow = c(1, 2))
  # plot(X)
  # plot(Y)
  
  # Basis expansions
  X_range <- X$rangeval[2] - X$rangeval[1]
  Y_range <- Y$rangeval[2] - Y$rangeval[1]
  
  basis_fourier_x <- t(eval.basis(evalarg = argvalsx,
                                  basisobj = create.fourier.basis(rangeval = c(X$rangeval[1], X$rangeval[2]), 
                                                                  nbasis = kx)))
  basis_fourier_y <- t(eval.basis(evalarg = argvalsy,
                                  basisobj = create.fourier.basis(rangeval = c(Y$rangeval[1], Y$rangeval[2]), 
                                                                  nbasis = ky)))
  
  for (i in 1:kx){
    basis_fourier_x[i,] <-
      basis_fourier_x[i,]/sqrt(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
  }
  
  for (i in 1:ky){
    basis_fourier_y[i,] <-
      basis_fourier_y[i,]/sqrt(integrateSimp1D(basis_fourier_y[i,]^2, c(Y_range)))
  }
  
  basis_fourier_x <- t(basis_fourier_x); basis_fourier_y <- t(basis_fourier_y)
  
  projX <- matrix(0, nrow = n, ncol = kx)
  projY <- matrix(0, nrow = n, ncol = ky)
  
  for (i in 1:n){
    for (j in 1:kx){
      projX[i,j] <- integrateSimp1D(X$data[i,]*basis_fourier_x[,j],X_range)
    }
  }
  
  for (i in 1:n){
    for (j in 1:ky){
      projY[i,j] <- integrateSimp1D(Y$data[i,]*basis_fourier_y[,j],Y_range)
    }
  }
  
  # Adot
  Adot.vec=Adot(X)

  # Obtain the entire matrix Adot
  Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
  Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
  Ad=t(Ad)
  Ad=Ad+t(Ad)-diag(diag(Ad))
  
  t[[counter]] <- test(projX, projY, B_boot <-  500, Ad)
  setTxtProgressBar(pb, counter / M)
  
}

# Density statistic vs densities bootstraped statistics
orig_stats <- sapply(t, function(x) x$orig_stat)
plot(density(orig_stats), lwd = 5, xlim = c(0, max(orig_stats)))
rug(orig_stats)
sapply(t[1:50], function(x) lines(density(x$boot_stat), 
                                  col = rgb(1, 0, 0, alpha = 0.25)))

# Distribution of p-values
hist(sapply(t, function(x) x$p_value))
