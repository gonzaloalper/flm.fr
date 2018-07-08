## FLM WITH FUNCTIONAL RESPONSE EXAMPLE

library(viridis)

# setwd("C:/Users/gonza/Desktop")
# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
lx <- 201
ly <- 101
ss <- seq(0, 50, l = lx) # Argvals of X
tt <- seq(1, 2, l = ly) # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  #5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
  cos(t * s / pi)^2
}

# Visualization
surface_beta <- outer(ss, tt, FUN = beta)

# Sample covariate
library(fda.usc)
n <- 200
X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)

# Error
noise <- matrix(rnorm(n = n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = noise, argvals = tt)
Y$data <- X$data %*% surface_beta + noise

# Plots
par(mfrow = c(1, 2))
plot(X)
plot(Y)

npcX<-200
npcY<-100
pcX <- fpc2(X,npcX,equispaced = TRUE)
pcY <- fpc2(Y,npcY,equispaced = TRUE)

hat_b <- t(pcX$x[,1:npcX])%*%pcX$x[,1:npcX]
hat_b <- pseudoinverse(hat_b)
hat_b <- hat_b%*%t(pcX$x[,1:npcX])
hat_b <- hat_b%*%pcY$x[,1:npcY]*sqrt(pcX$h)*sqrt(pcY$h)

print(hat_b)

#====================================================================

## SURFACES

surface_beta_hat <- matrix(0, nrow = lx, ncol = ly)
for (i in 1:npcX){
  for (j in 1:npcY){
    #acc<-hat_beta[i,j]*outer(pcX$rotation[,i], pcY$rotation[,j])*sqrt(pcX$h/pcY$h)
    acc <- hat_b[i,j]*pcX$rotation[,i]%*%t(pcY$rotation[,j])*sqrt(pcX$h/pcY$h)
    surface_beta_hat <- surface_beta_hat+acc
    # surface_beta <- outer(ss, tt, FUN = beta)
  }
}

limites = c(min(surface_beta_hat), max(surface_beta_hat))
image(ss, tt, surface_beta_hat,col = viridis(20))
image(ss, tt, surface_beta, col = viridis(20), zlim = limites)

#persp(fdataobj$argvals, Y$argvals, surface_beta_hat, phi = 0, theta = 60)
#persp(ss, tt, surface_beta, phi = 0, theta = 60)

#====================================================================

theoretical_beta <- matrix(0, nrow = npcX, ncol = npcY)
acc_aux<-outer(ss, tt, FUN = beta)
for (i in 1:npcX){
  for (j in 1:npcY){
    acc<-acc_aux*outer(pcX$rotation[,i],pcY$rotation[,j])*sqrt(pcY$h/pcX$h)
    theoretical_beta[i,j]<-sum(acc)*(pcX$h)*(pcY$h)
    #theoretical_beta[i,j]<-sdetorus::integrateSimp2D(acc,c(1,1))
  }
}

diferencia <- (sum(((hat_b - theoretical_beta)/theoretical_beta)^2))^(1/(npcX*npcY)); 
diferencia

explained_variance_X <- sum(pcX$d[1:npcX])/sum(pcX$d); explained_variance_X
explained_variance_Y <- sum(pcY$d[1:npcY])/sum(pcY$d); explained_variance_Y

hat_y <- pcX$x%*%hat_b
residuals <- abs(hat_y - pcY$x)

PCvM_statistic(residuals)
