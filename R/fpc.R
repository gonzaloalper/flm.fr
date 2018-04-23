fpc = function (fdataobj, ncomp = 3) 
{
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
  scores <- matrix(0, ncol = J, nrow = n)
  l <- 1:ncomp

  scores[, 1:J_min] <- inprod.fdata(X_cen.fdata, vs)
  colnames(scores) <- paste("PC", 1:J, sep = "")

  out <- list(rotation = vs[l], x = scores, fdataobj.cen = X_cen.fdata,
              mean = x_mean, l = l, u = u[, l, drop = FALSE])
  class(out) = "fdata.comp"
  return(out)
}
