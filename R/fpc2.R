fpc2 = function (fdataobj, ncomp = 3, equispaced) 
{
  X <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  nam <- fdataobj[["names"]]
  
  # Center data. fda.usc's is slower
  x_mean <- colMeans(fdataobj$data)
  X_cen <- t(t(fdataobj$data) - x_mean)
  
  # Componentes
  l <- 1:ncomp
  
  # SVD way
  eigenres <- svd(X_cen)
  v <- eigenres$v
  d <- eigenres$d
  #if(missing(equispaced)){
  #  #comprobar si lo son
  #  sep <- diff(X$argvals)
  #  equi <- all(sep[1] == sep[-1])
  #} 
  
  #else if(equispaced == TRUE){
    h <- (rtt[2] - rtt[1])/(length(tt) - 1)
    v <- v/sqrt(h) # functional normalization
    # Scores al estilo de componentes principales
    scores1 <- h * X_cen %*% v[, l]
  #}
  
  #else if(equispaced == FALSE){
    #v <- v/sqrt(h) # functional normalization
    # Scores al estilo de componentes principales
    #scores1 <- h * X_cen %*% v[, l]
  #}
  
  out <- list(d = d, rotation = v[, l], x = scores1, mean = x_mean, l = l, h = h)
  class(out) = "fdata.comp"
  return(out)
}
