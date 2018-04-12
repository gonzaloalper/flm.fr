
library(viridis)

# Argvals
l <- 201
ss <- seq(0, 1, l = l) # Argvals of X
tt <- ss # Argvals of Y

# Theoretical beta
beta <- function(s, t) {
  5 / (1 + (s - 0.5)^2) * cos(2 * pi * t * s)
  #(s + t^2) / 100
}

# Visualization
image(ss, tt, outer(ss, tt, FUN = beta), col = viridis(20))

# Sample covatiate
library(fda.usc)
n <- 100
X <- r.ou(n = n, x0 = seq(-5, 5, l = n), alpha = 2)

# Error
noise <- matrix(rnorm(n = n * l, mean = 0, sd = 0.1), nrow = n, ncol = l)

# Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
Y <- fdata(mdata = ss, argvals = ss)
Y$data <- X$data %*% outer(ss, tt, beta) + noise

# Plots
par(mfrow = c(1, 2))
plot(X)
plot(Y)


