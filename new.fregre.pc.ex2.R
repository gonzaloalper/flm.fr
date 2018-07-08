
library(fda.usc)
library(viridis)

# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

## Sample from the model

beta <- function(s, t) {
  
  #5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
  cos(t * s^2 * pi)^2 / 100
  
}

sample_flm_fr <- function(n, argvalsX, argvalsY, beta, sigma = 0.5) {
  
  # Sample covariate, a Brownian motion
  X <- rproc2fdata(n = n, t = argvalsX, sigma = "brownian")
  
  # Continuous error
  ly <- length(argvalsY)
  noise <- r.ou(n = n, x0 = rep(0, n), alpha = 2, 
                t = argvalsY, sigma = sigma)
  
  # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
  Y <- noise
  surface_beta <- outer(argvalsX, argvalsY, FUN = beta)
  Y$data <- X$data %*% surface_beta + noise$data
  
  # Return all useful objects
  return(list("X" = X, "Y" = Y, "error" = noise, "beta" = beta, 
              "surface_beta" = surface_beta))
  
}

# Test of the model
set.seed(1234567)
n <- 20
argvalsX <- seq(0, 1, l = 201)
argvalsY <- seq(0, 2, l = 301)
samp <- sample_flm_fr(n = n, argvalsX = argvalsX, argvalsY = argvalsY, 
                      beta = beta)
par(mfrow = c(2, 2))
col <- rainbow(n)[rank(samp$X$data[, 1])]
plot(samp$X, main = "X", col = col, lty = 1)
plot(samp$Y, main = "Y", col = col, lty = 1)
plot(samp$error, main = "error", col = col, lty = 1)
image(samp$X$argvals, samp$Y$argvals, samp$surface_beta, 
      col = viridis(20), main = "beta")

## Theoretical eigenfunctions of the covariance operator in [0, 1]
## (theoretical principal components)

eigen_brownian <- function(argvals, j) {
  # Check page 35 of Cuesta-Albertos et al. (2018, arXiv version)
  fdata(mdata = sqrt(2) * sin(((j - 0.5) * pi) %o% argvals), argvals = argvals)
}

set.seed(1234567)
large_samp <- sample_flm_fr(n = 1e4, argvalsX = argvalsX, argvalsY = argvalsY,
                            beta = beta)
large_PCX <- fpc2(fdataobj = large_samp$X, ncomp = 10, equispaced = TRUE)
large_PCY <- fpc2(fdataobj = large_samp$Y, ncomp = 10, equispaced = TRUE)
large_eigfun_X <- fdata(mdata = t(large_PCX$rotation), argvals = large_samp$X$argvals)
large_eigfun_Y <- fdata(mdata = t(large_PCY$rotation), argvals = large_samp$Y$argvals)
par(mfrow = c(1, 2))
plot(large_eigfun_X[1:5])
plot(large_eigfun_Y[1:5])

# Correct estimation of the true theoretical eigenfunctions of the Brownian motion
par(mfrow = c(1, 3))
plot(eigen_brownian(argvals = argvalsX, j = 1:5), main = "Theoretical") 
plot(large_eigfun_X[1:5] * sign(large_eigfun_X$data[1:5, 2]), main = "Empirical") # Change signs
plot(eigen_brownian(argvals = argvalsX, j = 1:5), main = "Both")
lines(large_eigfun_X[1:5] * sign(large_eigfun_X$data[1:5, 2]))

## Representation of data in the eigenfunctions (both theoretical and estimated)

sqrt(diff(argvalsX)[1] * rowSums(eigen_brownian(argvals = argvalsX, j = 1:5)$data^2)) # Unit functional norm


s <- (sign(eigen_brownian(argvals = argvalsX, j = 1:5)[, 10]$data) != 
        sign(large_eigfun_X[1:5]$data[, 10]))
s[s == 1] <- -1
s[s == 0] <- 1

# Theoretical scores
theoScores <- matrix(nrow = 20, ncol = 5)
for (i in 1:20) {
  for (j in 1:5) {
    theoScores[i, j] <- diff(argvalsX)[1] * rowSums(eigen_brownian(argvals = argvalsX, j = j)$data * large_samp$X$data[i, ]) # Unit functional norm
  }
}

# Empirical scores
empScores <- fpc2(fdataobj = large_samp$X, ncomp = 5)$x[1:20, ]
dif <- theoScores - empScores * matrix(s, byrow = TRUE, nrow = 20, ncol = 5)
max(abs(dif))

# fpc2 works as expected

## Some elements of the tensor basis

par(mfrow = c(3, 3))
for (i in 1:3) {
  for (j in 1:3) {
    image(large_samp$X$argvals, large_samp$Y$argvals, 
          outer(large_eigfun_X$data[i, ], large_eigfun_Y$data[j, ]), 
          col = viridis(20), main = paste0("PC", i, " x ", "PC", j))
  }
}

## Projection of beta to the (kx, ky)-truncated tensor basis

beta_proj_basis <- function(surface_beta, eigfun_X, eigfun_Y, kx, ky) {
  
  surface_beta_k <- 0
  B <- matrix(nrow = kx, ncol = ky)
  hx <- diff(eigfun_X$argvals[1:2])
  hy <- diff(eigfun_Y$argvals[1:2])
  for (i in 1:kx) {
    for (j in 1:ky) {
      
      PCij <- outer(eigfun_X$data[i, ], eigfun_Y$data[j, ])
      B[i, j] <- sum(surface_beta * PCij) * hx * hy
      surface_beta_k <- surface_beta_k + B[i, j] * PCij
      
    }
    
  }
  
  return(list("B" = B, "surface_beta" = surface_beta_k))
  
}

proj1 <- beta_proj_basis(surface_beta = large_samp$surface_beta, 
                         eigfun_X = large_eigfun_X, eigfun_Y = large_eigfun_Y, 
                         kx = 1, ky = 2)
proj2 <- beta_proj_basis(surface_beta = large_samp$surface_beta, 
                         eigfun_X = large_eigfun_X, eigfun_Y = large_eigfun_Y, 
                         kx = 2, ky = 3)
proj3 <- beta_proj_basis(surface_beta = large_samp$surface_beta, 
                         eigfun_X = large_eigfun_X, eigfun_Y = large_eigfun_Y, 
                         kx = 3, ky = 3)
proj5 <- beta_proj_basis(surface_beta = large_samp$surface_beta, 
                         eigfun_X = large_eigfun_X, eigfun_Y = large_eigfun_Y, 
                         kx = 5, ky = 5)
proj10 <- beta_proj_basis(surface_beta = large_samp$surface_beta, 
                          eigfun_X = large_eigfun_X, eigfun_Y = large_eigfun_Y, 
                          kx = 10, ky = 10)
par(mfrow = c(3, 2))
lim <- c(-0.005, 0.0175)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, proj1$surface_beta, 
      col = viridis(20), main = "beta_{1, 2}", zlim = lim)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, proj2$surface_beta, 
      col = viridis(20), main = "beta_{2, 3}", zlim = lim)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, proj3$surface_beta, 
      col = viridis(20), main = "beta_{3, 3}", zlim = lim)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, proj5$surface_beta, 
      col = viridis(20), main = "beta_{5, 5}", zlim = lim)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, proj10$surface_beta, 
      col = viridis(20), main = "beta_{10, 10}", zlim = lim)
image(large_eigfun_X$argvals, large_eigfun_Y$argvals, large_samp$surface_beta, 
      col = viridis(20), main = "beta_{inf, inf}", zlim = lim)

# Increasing norm explanation
sum(proj1$surface_beta^2)
sum(proj2$surface_beta^2)
sum(proj3$surface_beta^2)
sum(proj5$surface_beta^2)
sum(proj10$surface_beta^2)
sum(large_samp$surface_beta^2)

## Estimation of the (kx, ky)-projected beta

set.seed(1234567)
n <- 200
small_samp <- sample_flm_fr(n = n, argvalsX = seq(0, 1, l = 201), 
                            argvalsY = seq(0, 2, l = 301), beta = beta)
kx <- 5
ky <- 5
small_PCX <- fpc2(fdataobj = small_samp$X, ncomp = kx, equispaced = TRUE)
small_PCY <- fpc2(fdataobj = small_samp$Y, ncomp = ky, equispaced = TRUE)
small_eigfun_X <- fdata(mdata = t(small_PCX$rotation), argvals = small_samp$X$argvals)
small_eigfun_Y <- fdata(mdata = t(small_PCY$rotation), argvals = small_samp$Y$argvals)
par(mfrow = c(2, 2))
plot(large_eigfun_X[1:5], ylim = c(-2, 2), main = "X theoretical")
plot(small_eigfun_X[1:5], ylim = c(-2, 2), main = "X estimated")
plot(large_eigfun_Y[1:5], ylim = c(-2, 2), main = "Y theoretical")
plot(small_eigfun_Y[1:5], ylim = c(-2, 2), main = "Y estimated")
# Caution with the signs!

# NO OK
pseudoinverse(crossprod(small_PCX$x)) %*% t(small_PCX$x) %*% small_PCY$x
proj5$B

# OK
pseudoinverse(crossprod(small_PCX$x)) %*% t(small_PCX$x) %*% small_PCY$x * 
  sqrt(small_PCX$h / small_PCY$h)
proj5$B / sqrt(small_PCX$h * small_PCY$h)
# Factor mismatch. Reason must be in the linear model formulation!
