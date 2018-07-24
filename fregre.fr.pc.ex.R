## FLM WITH FUNCTIONAL RESPONSE EXAMPLE

library(viridis)

setwd("C:/Users/gonza/Desktop")
# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
lx <- 101
ly <- 101
ss <- seq(3, 4, l = lx) # Argvals of X
tt <- seq(1, 2, l = ly) # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  #diag(2,lx,ly)
  #2+(s - s + t - t)
  #5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  (s + t^2) / 100
  #cos(t * s / pi)^2
}

# Visualization
surface_beta <- outer(ss, tt, FUN = beta)

# Sample covariate
library(fda.usc)
n <- 200
X <- rproc2fdata(n = n, ss, sigma="brownian")
#X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)
X_range <- X$rangeval[2] - X$rangeval[1]

# Error
#noise <- matrix(0, nrow = n, ncol = ly)
noise <- matrix(rnorm(n = n * ly, mean = 0, sd = 0), nrow = n, ncol = ly)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = noise, argvals = tt)

#Y$data <- X$data %*% surface_beta * hx

for (j in 1:ly){
  b = surface_beta[,j]
  sapply(1:n, function(i){
    Y$data[i,j] <<- sdetorus::integrateSimp1D(b * X$data[i,], X_range)
  })
}

Y$data <- Y$data + noise 

Y_range <- Y$rangeval[2] - Y$rangeval[1]

# Plots
par(mfrow = c(1, 2))
plot(X)
plot(Y)

npcX<-10
npcY<-10
pcX <- fpc2(X,npcX,equispaced = TRUE)
pcY <- fpc2(Y,npcY,equispaced = TRUE)

hat_b <- t(pcX$x[,1:npcX])%*%pcX$x[,1:npcX]
hat_b <- pseudoinverse(hat_b)
hat_b <- hat_b%*%t(pcX$x[,1:npcX])
hat_b <- hat_b%*%pcY$x[,1:npcY]

print(hat_b)

#====================================================================

## SURFACES

surface_beta_hat <- matrix(0, nrow = lx, ncol = ly)
for (i in 1:npcX){
  for (j in 1:npcY){
    acc<-hat_b[i,j]*outer(pcX$rotation[,i], pcY$rotation[,j])
    #acc <- hat_b[i,j]*pcX$rotation[,i]%*%t(pcY$rotation[,j])
    surface_beta_hat <- surface_beta_hat+acc
    # surface_beta <- outer(ss, tt, FUN = beta)
  }
}

minima <- c(min(surface_beta_hat), min(surface_beta))
maxima <- c(max(surface_beta_hat), max(surface_beta))
limits = c(min(minima), max(maxima))
image(ss, tt, surface_beta_hat, col = viridis(20), zlim = limits)
image(ss, tt, surface_beta, col = viridis(20), zlim = limits)

#persp(fdataobj$argvals, Y$argvals, surface_beta_hat, phi = 0, theta = 60)
#persp(ss, tt, surface_beta, phi = 0, theta = 60)

#====================================================================

theoretical_beta <- matrix(0, nrow = npcX, ncol = npcY)
acc_aux<-outer(ss, tt, FUN = beta)
for (i in 1:npcX){
  for (j in 1:npcY){
    acc<-acc_aux*outer(pcX$rotation[,i],pcY$rotation[,j])
    theoretical_beta[i,j]<-sum(acc)*pcX$h*pcY$h
    #theoretical_beta[i,j]<-sdetorus::integrateSimp2D(acc,c(1,1))
  }
}

t_beta <- matrix(0, nrow = npcX, ncol = npcY)
for (i in 1:npcX){
  for (j in 1:npcY){
    feval <- acc_aux[i,j] * 2 * outer(sin((i - 0.5) * pi * ss),sin((j - 0.5) * pi * tt))
    t_beta[i,i] <- sdetorus::integrateSimp2D(feval,c(1,1))
  }
}

diferencia <- (sum(((hat_b - theoretical_beta)/theoretical_beta)^2))^(1/(npcX*npcY)); 
diferencia

explained_variance_X <- sum(pcX$d[1:npcX])/sum(pcX$d); explained_variance_X
explained_variance_Y <- sum(pcY$d[1:npcY])/sum(pcY$d); explained_variance_Y

pvalue <- flm_test(10000, n, X, pcX, pcY, hat_b)
