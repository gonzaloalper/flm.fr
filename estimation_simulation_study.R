library(fda.usc)
library(viridis)
#library(devtools)
#install_github("egarpor/rp.flm.test")
library(rp.flm.test)
# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
lx <- 201
ly <- 201
ss <- seq(0, 2, l = lx) # Argvals of X
tt <- seq(0, 3, l = ly) # Argvals of Y

# Sample covariate
n <- 80
X <- rproc2fdata(n = n, ss, sigma="brownian")
#X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)
X_range <- X$rangeval[2] - X$rangeval[1]
hx <- (X$rangeval[2]-X$rangeval[1])/length(ss)

# Theoretical beta
beta <- function(s, t) {
  #diag(2,lx,ly)
  #2+(s - s + t - t)
  #cos(2 * pi * t * s) / (1 + (s - 0.5)^2) 
  #1 - ((s - 2) - (t - 3))^2
  sin(6*pi*s) + cos(6*pi*t)
  #pi^2*(s^2-1/(0.25+t)) * tanh(1+s + t^3)
  #cos(t * s / pi)^2
}

# Visualization
surface_beta <- outer(ss, tt, FUN = beta)

Y <- linear_model(X, tt, 0.01, beta)

Y_range <- Y$rangeval[2] - Y$rangeval[1]

# Center data. fda.usc's is slower
X_cen <- X
x_mean <- colMeans(X$data)
X_cen$data <- t(t(X$data) - x_mean)

Y_cen <- Y
y_mean <- colMeans(Y$data)
Y_cen$data <- t(t(Y$data) - y_mean)

# Plots
par(mfrow = c(2, 2))
plot(X)
plot(Y)
plot(X_cen)
plot(Y_cen)

npcX <- 73
npcY <- 2
pcX <- fpc(X_cen,npcX,equispaced = TRUE)
pcY <- fpc(Y_cen,npcY,equispaced = TRUE)

# RSSN:
hat_b_aux <- t(pcX$x[,1:npcX])%*%pcX$x[,1:npcX]
hat_b_aux <- pseudoinverse(hat_b_aux)
hat_b_aux <- hat_b_aux %*% t(pcX$x[,1:npcX])
hat_b <- hat_b_aux %*% pcY$x[,1:npcY]

# Regularization:
#hat_b <- beta_hat_net(pcX$x, pcY$x, type = "lasso", k_folds = n)$beta_hat

theoretical_beta <- matrix(0, nrow = npcX, ncol = npcY)

for (i in 1:npcX){
  for (j in 1:npcY){
    outer_product <- outer(pcX$rotation[,i],pcY$rotation[,j])
    acc <- surface_beta * outer_product
    theoretical_beta[i,j] <- integrateSimp2D(acc, c(X_range, Y_range))
  }
}

#EXPLAINED VARIANCE

diferencia <- (sum(((hat_b - theoretical_beta)/theoretical_beta)^2))^(1/(npcX*npcY)); 
explained_variance_X <- sum(pcX$d[1:npcX])/sum(pcX$d); explained_variance_X
explained_variance_Y <- sum(pcY$d[1:npcY])/sum(pcY$d); explained_variance_Y

## SURFACES
surface_beta_hat <- matrix(0, nrow = lx, ncol = ly)
surface_pca <- matrix(0, nrow = lx, ncol = ly)
for (i in 1:npcX){
  for (j in 1:npcY){
    outer_product <- outer(pcX$rotation[,i], pcY$rotation[,j])
    acc1 <- theoretical_beta[i,j] * outer_product
    surface_pca <- surface_pca + acc1
    acc2 <- hat_b[i,j] * outer_product
    surface_beta_hat <- surface_beta_hat + acc2
  }
}

par(mfrow = c(1, 3))
minima <- c(min(surface_beta), min(surface_pca), min(surface_beta_hat))
maxima <- c(max(surface_beta), max(surface_pca), max(surface_beta_hat))
limits = c(min(minima), max(maxima))
image(ss, tt, surface_beta, col = viridis(20), zlim = limits, xlab = "X", ylab = "Y")
image(ss, tt, surface_pca, col = viridis(20), zlim = limits, xlab = "X", ylab = "Y")
image(ss, tt, surface_beta_hat, col = viridis(20), zlim = limits, xlab = "X", ylab = "Y")

acc <- 0
for (i in 1:201){
  for (j in 1:201){
    acc <- acc + (surface_beta_hat[i,j]-surface_beta[i,j])^2
  }
}
sqrt(acc/201^2)