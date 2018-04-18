 # Functional Model with Scalar Response
 # Regression on PCs

fregre.sr.pc=function (fdataobj, Y, ncomp = 3) 
{
  #==============================
  # COMPUTE PRINCIPAL COMPONENTS
  #==============================
  
  C <- match.call()
  if (!is.fdata(fdataobj)) 
    stop("No fdata class")
  nas1 <- is.na.fdata(fdataobj)
  if (any(nas1)) 
    stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
  X <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  nam <- fdataobj[["names"]]
  mm <- fdata.cen(fdataobj)
  xmean <- mm$meanX
  Xcen.fdata <- mm$Xcen
  dimx <- dim(X)
  n <- dimx[1]
  J <- dimx[2]
  Jmin <- min(c(J, n))
  M <- solve(diag(J))
  Xcen.fdata$data <- Xcen.fdata$data %*% t(M)
  eigenres <- svd(Xcen.fdata$data)
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
  scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs)
  
  colnames(scores) <- paste("PC", 1:J, sep = "")
  l <- 1:ncomp
  
  #==============================
  # REGRESSION
  #==============================
  
  object.lm<-lm(Y~scores[,l])
  
  object.lm = lm(Y~scores[,l])
  e<-object.lm$residuals
  
  out <- list(call = C, d = newd, rotation = vs[1:ncomp], x = scores, 
              fdataobj.cen = Xcen.fdata, norm = norm, type = "pc", mean = xmean, 
              fdataobj = fdataobj, l = l, u = u[, 1:ncomp, drop = FALSE],
              coefficients=object.lm$coefficients,residuals = e,
              fitted.values =object.lm$fitted.values)
  class(out) = "fdata.comp"
  return(out)
}
