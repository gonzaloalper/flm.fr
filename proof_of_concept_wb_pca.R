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
X <- r.ou(n = n, x0 = seq(-1, 1, l = n), alpha = 1/3, t = argvalsx)

# Theoretical beta
beta <- function(s, t) {
  pi^2*(s^2-1/(0.25+t)) * tanh(1+s + t^3)
}

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
d <- 0.01
Y1 <- Y
Y$data <- Y$data + noise + d * exp(X$data)
Y1$data <- Y$data + exp(X$data) 

# # Plots
# X$data <- X$data[1:10,]
# Y$data <- Y$data[1:10,]
par(mfrow = c(1, 2))
plot(X,xlab = "s", ylab = "X(s)",ylim=c(-5,5))
plot(Y1,xlab = "t", ylab = "Y(t)",ylim=c(-100,100))

# Basis expansions
kx <- 5; ky <- 5
X_range <- X$rangeval[2] - X$rangeval[1]
Y_range <- Y$rangeval[2] - Y$rangeval[1]

pcX <- fpc(X,kx,equispaced = TRUE)
pcY <- fpc(Y,ky,equispaced = TRUE)

projX <- pcX$x; projY <- pcY$x

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
pb <- txtProgressBar(style = 3)kx

for (counter in 1:M) {
  # Data generation
  n <- 200
  kx <- 5; ky <- 5
  # Data generation
  X <- r.ou(n = n, x0 = seq(-1, 1, l = n), alpha = 1/3, t = argvalsx)
  
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
  d <- 0.010
  Y$data <- Y$data + noise + d * exp(X$data)
  
  # # Plots
  # par(mfrow = c(1, 2))
  # plot(X)
  # plot(Y)
  
  # Basis expansions
  X_range <- X$rangeval[2] - X$rangeval[1]
  Y_range <- Y$rangeval[2] - Y$rangeval[1]
  
  pcX <- fpc(X,kx,equispaced = TRUE)
  pcY <- fpc(Y,ky,equispaced = TRUE)
  
  projX <- pcX$x; projY <- pcY$x
  
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
