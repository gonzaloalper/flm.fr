fpc3 = function (fdataobj, ncomp = 3) 
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
  
  eigenres <- svd(X_cen)
  v <- eigenres$v

  # Scores usando inprod.fdata (aproxima el producto interior como se hace en fda.usc)
  vs <- fdata(t(v[, l]), argvals = tt, rtt, list(main = "fdata2pc", xlab = "t", 
                                               ylab = "rotation"))
  scores2 <- length(tt)/rtt[2]*inprod.fdata(fdataobj, vs) # Divide por el paso de integración!


  out <- list(rotation = v[, l], x = scores2, mean = x_mean, l = l)
  class(out) = "fdata.comp"
  return(out)
}
