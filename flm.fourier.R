## FLM WITH FUNCTIONAL RESPONSE EXAMPLE

library(viridis)
library(fda.usc)

setwd("C:/Users/gonza/Desktop")
# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
lx <- 201
ly <- 201
ss <- seq(3, 5, l = lx) # Argvals of X
tt <- seq(1, 2, l = ly) # Argvals of Y
n_basis_x <- 15; n_basis_y <- 15

# Sample covariate
n <- 2000
X <- rproc2fdata(n = n, ss, sigma="brownian")
#X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)
X_range <- X$rangeval[2] - X$rangeval[1]
hx <- (X$rangeval[2]-X$rangeval[1])/length(ss)

# Theoretical beta
beta <- function(s, t) {
  #diag(2,lx,ly)
  #2+(s - s + t - t)
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
  #cos(t * s / pi)^2
}

# Visualization
surface_beta <- outer(ss, tt, FUN = beta)

Y <- linear_model(X, n, ly, tt, 0, surface_beta)

Y_range <- Y$rangeval[2] - Y$rangeval[1]

# Plots
par(mfrow = c(1, 2))
plot(X)
plot(Y)

basis_fourier_x <- t(eval.basis(evalarg = seq(X$rangeval[1], X$rangeval[2], l = lx),
                                basisobj = create.fourier.basis(rangeval = c(X$rangeval[1], X$rangeval[2]), 
                                                                nbasis = n_basis_x)))
basis_fourier_y <- t(eval.basis(evalarg = seq(Y$rangeval[1], Y$rangeval[2], l = ly),
                                basisobj = create.fourier.basis(rangeval = c(Y$rangeval[1], Y$rangeval[2]), 
                                                                nbasis = n_basis_y)))

for (i in 1:n_basis){
  print(sdetorus::integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
}

for (i in 1:n_basis_x){
  basis_fourier_x[i,] <-
    basis_fourier_x[i,]/sqrt(sdetorus::integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
}

for (i in 1:n_basis_x){
  print(sdetorus::integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
}

for (i in 1:n_basis_y){
  basis_fourier_y[i,] <-
    basis_fourier_y[i,]/sqrt(sdetorus::integrateSimp1D(basis_fourier_y[i,]^2, c(Y_range)))
}

basis_fourier_x <- t(basis_fourier_x); basis_fourier_y <- t(basis_fourier_y)

par(mfrow = c(1, 1)); matplot(basis_fourier_x, type="l")

projX <- matrix(0, nrow = n, ncol = n_basis_x)
projY <- matrix(0, nrow = n, ncol = n_basis_y)

for (i in 1:n){
  for (j in 1:n_basis_x){
    projX[i,j] <- sdetorus::integrateSimp1D(X$data[i,]*basis_fourier_x[,j],X_range)
  }
}

for (i in 1:n){
  for (j in 1:n_basis_y){
    projY[i,j] <- sdetorus::integrateSimp1D(Y$data[i,]*basis_fourier_y[,j],Y_range)
  }
}

hat_b <- pseudoinverse(t(projX) %*% projX) %*% t(projX) %*% projY

coefs_beta <- matrix(0, nrow = n_basis_x, ncol = n_basis_y)

for (i in 1:n_basis_x){
  for (j in 1:n_basis_y){
    outer_product <- outer(basis_fourier_x[,i],basis_fourier_y[,j])
    acc <- surface_beta * outer_product
    coefs_beta[i,j] <- sdetorus::integrateSimp2D(acc, c(X_range, Y_range))
  }
}

sdetorus::integrateSimp2D(outer(basis_fourier_x[,i],basis_fourier_y[,j])^2,
                          c(X_range,Y_range))

surface_fourier_hat <- matrix(0, nrow = lx, ncol = ly)
surface_fourier <- matrix(0, nrow = lx, ncol = ly)

for (i in 1:n_basis_x){
  for (j in 1:n_basis_y){
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
image(ss, tt, surface_fourier_hat, col = viridis(20), zlim = limits)
image(ss, tt, surface_fourier, col = viridis(20), zlim = limits)
image(ss, tt, surface_beta, col = viridis(20), zlim = limits)
