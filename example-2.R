#Datos ejemplo

# Argvals
lx <- 201
ly <- 101
ss <- seq(0, 1, l = lx) # Argvals of X
tt <- seq(0, 1, l = ly) # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
}

# Visualization
image(ss, tt, outer(ss, tt, FUN = beta), col = viridisLite::viridis(20))

# Sample covariate
library(fda.usc)
n <- 20
fdataobj <- r.ou(n = n, x0 = seq(-5, 5, l = n), alpha = 2, t = ss)

# Error
noise <- matrix(rnorm(n = n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = noise, argvals = tt)
Y$data <- fdataobj$data %*% outer(ss, tt, beta) + noise

#==========================================================

# PRINCIPAL COMPONENTS

ncomp <- 3

#browser()
X <- fdataobj[["data"]]
tt <- fdataobj[["argvals"]]
rtt <- fdataobj[["rangeval"]]
nam <- fdataobj[["names"]]

# Center data. fda.usc's is slower
x_mean <- colMeans(fdataobj$data)
X_cen <- t(t(fdataobj$data) - x_mean)
mm <- fdata.cen(fdataobj)

# Necesario?  
# dim_x <- dim(X)
# n <- dim_x[1]
# J <- dim_x[2]
# J_min <- min(c(J, n))

# Componentes
ncomp <- 3
l <- 1:ncomp

# These are the eigenfunctions! They are vectors of unit norm

# SVD way
eigenres <- svd(X_cen)
v <- eigenres$v
# u <- eigenres$u # No necesario?
d <- eigenres$d

# Eigendecomposition of the covariance matrix way
eig <- eigen(crossprod(X_cen))

# Compare
d^2
eig$values

# Same orthogonal vectors (of unit norm)
eig$vectors[, l]
v[, l]
colSums(eig$vectors[, l]^2)

# Scores al estilo de componentes principales
scores1 <- X %*% v[, l]
scores1

# Scores usando inprod.fdata (aproxima el producto interior como se hace en fda.usc)
vs <- fdata(t(v[, l]), argvals = tt, rtt, list(main = "fdata2pc", xlab = "t", 
                                          ylab = "rotation"))
scores2 <- inprod.fdata(fdataobj, vs) # Divide por el paso de integración!
scores2
# Aprox integral = sum / paso integracion, por eso la diferencia en escala
# La diferencia en escala no es exacta porque se usa una cuadratura de Simpson y 
# ésta otorga pesos distintos a la discretización (particularmente, en los extremos)
# Pero con lx elevado es bastante similar a 1 / (length(ss) - 1)
colMeans(abs(scores2 / scores1))
1 / (lx - 1)

# Estos dos tipos de scores corresponden a dos de las 3 variantes descritas en el
# libro de Ramsay y Silverman. La primera (pasar al caso multivariante sin más)
# sería adecuada sólamente si el grid fuese equiespaciado, pero hay que mirarlo. 
# La segunda utiliza integración numérica para aproximar el producto interior de
# la proyección

# Comprobación de la aproximación de los datos por ambos scores
par(mfrow = c(1, 2))
matplot(t(X), type = "l", lty = 3, col = 1)
matlines(t(scores1 %*% t(v[, l])), type = "l", col = 2)
matplot(t(X), type = "l", lty = 3, col = 1)
matlines(t(scores2 %*% t(v[, l])), type = "l", col = 2) # Necesita re-escalar las v's, por eso la modificación en fda.usc!

