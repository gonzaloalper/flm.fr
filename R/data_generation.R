data_generation = function(n, ss, tt)
{
  # Sample covariate
  X <- rproc2fdata(n = n, ss, sigma="brownian")
  #X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)

  # Theoretical beta
  beta <- function(s, t) {
    #diag(2,lx,ly)
    #(s - s + t - t)
    #pi^2*(s^2-1/(0.25+t)) * tanh(1+s + t^2)
    #1 + sin(2*pi*s) + cos(2*pi*t) + sin(2*pi*s)*cos(2*pi*t)
    cos(2 * pi * t * s) / (1 + (s - 0.5)^2)
    #(s + t^2) / 100
    #cos(t * s / pi)^2
  }
  
  # Visualization
  surface_beta <- outer(ss, tt, FUN = beta)
  Y <- linear_model(X, tt, 0.01, beta)
  
  out <- list(X, Y, beta, surface_beta)
  
  return(out)
}
