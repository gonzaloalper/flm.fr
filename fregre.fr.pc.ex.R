
## FLM WITH FUNCTIONAL RESPONSE EXAMPLE

library(viridis)

# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
lx <- 201
ly <- 101
ss <- seq(0, 1, l = lx) # Argvals of X
tt <- seq(0, 2, l = ly) # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
}

# Visualization
image(ss, tt, outer(ss, tt, FUN = beta), col = viridis(20))

# Sample covariate
library(fda.usc)
n <- 200
fdataobj <- r.ou(n = n, x0 = seq(-5, 5, l = n), alpha = 2, t = ss)

# Error
noise <- matrix(rnorm(n = n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = noise, argvals = tt)
Y$data <- fdataobj$data %*% outer(ss, tt, beta) + noise

# Plots
par(mfrow = c(1, 2))
plot(fdataobj)
plot(Y)

npcX<-5
npcY<-5
pcX <- fpc2(fdataobj,npcX)
pcY <- fpc3(Y,npcY)

hat_beta <- t(pcX$x[,1:npcX])%*%pcX$x[,1:npcX]
hat_beta <- pseudoinverse(hat_beta)
hat_beta <- hat_beta%*%t(pcX$x[,1:npcX])
hat_beta <- hat_beta%*%pcY$x[,1:npcY]

print(hat_beta)

surface_beta_hat <- 0
for (i in 1:npcX){
  for (j in 1:npcY){
    acc<-hat_beta[i,j]*outer(pcX$rotation[,i],pcY$rotation[,j])
    surface_beta_hat<-surface_beta_hat+acc
  }
}

image(fdataobj$argvals, Y$argvals, surface_beta_hat,col = viridis(20))
image(ss, tt, surface_beta<-outer(ss, tt, FUN = beta), col = viridis(20))

theoretical_beta <- matrix(0, nrow = npcX, ncol = npcY)
acc_aux<-outer(ss, tt, FUN = beta)
for (i in 1:npcX){
  for (j in 1:npcY){
    acc<-acc_aux*outer(pcX$rotation[,i],pcY$rotation[,j])
    theoretical_beta[i,j]<-sum(acc)
    #theoretical_beta[i,j]<-sdetorus::integrateSimp2D(acc,c(1,1))
  }
}

abs(theoretical_beta-hat_beta)/theoretical_beta

theoretical_beta
hat_beta
