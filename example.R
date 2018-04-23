#Datos ejemplo

# Argvals
lx <- 11
ly <- 101
ss <- seq(0, 1, l = lx) # Argvals of X
tt <- seq(0, 1, l = ly) # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s) / 200
  #(s + t^2) / 100
}

# Visualization
image(ss, tt, outer(ss, tt, FUN = beta), col = viridis(20))

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
mm <- fdata.cen(fdataobj)

x_mean <- mm$meanX
X_cen.fdata <- mm$Xcen
dim_x <- dim(X)
n <- dim_x[1]
J <- dim_x[2]
J_min <- min(c(J, n))

eigenres <- svd(X_cen.fdata$data)
v <- eigenres$v
u <- eigenres$u
d <- eigenres$d
D <- diag(d)
vs <- fdata(t(v), argvals = tt, rtt, list(main = "fdata2pc", xlab = "t", 
                                          ylab = "rotation"))
scores <- matrix(0, ncol = ncomp, nrow = n)
l <- 1:ncomp

#scores[, 1:J_min] <- inprod.fdata(X_cen.fdata, vs[l,])

auto <- vs[l,]
mdist = array(0, dim = c(numgr, numgr2))

numgr <- nrow(X_cen.fdata)
numgr2 <- nrow(auto)

for (i in 1:numgr) {
  for (ii in 1:numgr2) {
    f = X_cen.fdata[i, ] * auto[ii, ]
    scores[i, ii] = fda.usc:::int.simpson2(tt, f$data)
  }
}

scores
#plot(fdataobj)
plot(fdata(scores[,l]%*%vs$data[l,]))
