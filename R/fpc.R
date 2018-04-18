fpc = function (fdataobj, ncomp = 3) 
{
  # Computes the principal components
  C <- match.call()
  if (!is.fdata(fdataobj)) 
    stop("No fdata class")
  nas1 <- is.na.fdata(fdataobj)
  if (any(nas1)) 
    stop("fdataobj contain ", sum(nas1), " There are curves with some NA value \n")

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

  M <- solve(diag(J))
  X_cen.fdata$data <- X_cen.fdata$data %*% t(M)
  eigenres <- svd(X_cen.fdata$data)
  v <- eigenres$v
  u <- eigenres$u
  d <- eigenres$d
  D <- diag(d)
  vs <- fdata(t(v), tt, rtt, list(main = "fdata2pc", xlab = "t", 
                                ylab = "rotation"))
  scores <- matrix(0, ncol = J, nrow = n)

  #Normalize the data:
  dtt <- diff(tt)
  drtt <- diff(rtt)
  eps <- as.double(.Machine[[1]] * 10)
  inf <- dtt - eps
  sup <- dtt + eps

  if (all(dtt > inf) & all(dtt < sup)) 
    delta <- sqrt(drtt/(J - 1))
  else delta <- 1/sqrt(mean(1/dtt))
  no <- norm.fdata(vs)
  vs <- vs/delta
  newd <- d * delta
  scores[, 1:J_min] <- inprod.fdata(X_cen.fdata, vs)

  colnames(scores) <- paste("PC", 1:J, sep = "")
  l <- 1:ncomp

  out <- list(rotation = vs[l], x = scores, fdataobj.cen = X_cen.fdata,
              mean = x_mean, l = l, u = u[, l, drop = FALSE])
  class(out) = "fdata.comp"
  return(out)
}
