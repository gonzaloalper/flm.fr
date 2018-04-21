
## FLM WITH FUNCTIONAL RESPONSE EXAMPLE

library(viridis)

# Execute all the functions in the R directory from the above file:
sourceDir <- function(directory = "R", ...) {
  invisible(sapply(dir(directory),
                   function(x) source(paste0(directory, "/", x), ...)))
}
sourceDir()

# Argvals
l <- 201
ss <- seq(0, 1, l = l) # Argvals of X
tt <- ss # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
}

# Visualization
image(ss, tt, outer(ss, tt, FUN = beta), col = viridis(20))

# Sample covatiate
library(fda.usc)
n <- 1000
fdataobj <- r.ou(n = n, x0 = seq(-5, 5, l = n), alpha = 2)

# Error
noise <- matrix(rnorm(n = n * l, mean = 0, sd = 0.1), nrow = n, ncol = l)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = ss, argvals = ss)
Y$data <- fdataobj$data %*% outer(ss, tt, beta) + noise

# Plots
par(mfrow = c(1, 2))
plot(fdataobj)
plot(Y)

npc<-5
pcX <- fpc(fdataobj,npc)
pcY <- fpc(Y,npc)

hat_beta <- t(pcX$x[,1:npc])%*%pcX$x[,1:npc]
hat_beta <- fda.usc:::Minverse(hat_beta)
hat_beta <- hat_beta%*%t(pcX$x[,1:npc])
hat_beta <- hat_beta%*%pcY$x[,1:npc]

print(hat_beta)

surface_beta<-0
for (i in 1:npc){
  for (j in 1:npc){
    acc<-hat_beta[i,j]*outer(pcX$rotation$data[i,],pcY$rotation$data[j,])
    surface_beta<-surface_beta+acc
  }
}

image(surface_beta,col = viridis(20))
image(ss, tt, outer(ss, tt, FUN = beta), col = viridis(20))

theoretical_beta <- matrix(0, nrow = npc, ncol = npc)
acc_aux<-outer(ss, tt, FUN = beta)
for (i in 1:npc){
  for (j in 1:npc){
    acc<-acc_aux%*%outer(pcX$rotation$data[i,],pcY$rotation$data[j,])
    theoretical_beta[i,j]<-sdetorus::integrateSimp2D(acc,c(1,1))
  }
}
