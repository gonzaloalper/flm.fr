fpc2 = function (fdataobj, ncomp = 3) 
{
  #browser()
  X <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  nam <- fdataobj[["names"]]
  
  # Center data. fda.usc's is slower
  x_mean <- colMeans(fdataobj$data)
  X_cen <- t(t(fdataobj$data) - x_mean)
  mm <- fdata.cen(fdataobj)
  
  # Componentes
  l <- 1:ncomp
  
  # SVD way
  eigenres <- svd(X_cen)
  v <- eigenres$v
  # u <- eigenres$u # No necesario?
  d <- eigenres$d
  
  # Eigendecomposition of the covariance matrix way
  eig <- eigen(crossprod(X_cen))
  
  # Scores al estilo de componentes principales
  scores1 <- X %*% v[, l]
  
  out <- list(rotation = v[, l], x = scores1, mean = x_mean, l = l)
  class(out) = "fdata.comp"
  return(out)
}
