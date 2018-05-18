PCvM_statistic = function(residuals)
  {
  library(fda.usc)
  
  # Adot
  Adot.vec=Adot(X)
  
  # Obtain the entire matrix Adot
  Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
  Ad[upper.tri(Ad,diag = FALSE)]=Adot.vec[-1]
  Ad=t(Ad)
  Ad=Ad+t(Ad)-diag(diag(Ad))
  
  n <- dim(residuals)[1]
  ky <- dim(residuals)[2]
  
  PCvM <- 0
  for (i in 1:ky){
    acc <- t(residuals)[i,]%*%Ad%*%residuals[,i]
    PCvM <- PCvM + acc
  }
  
  PCvM <- pi^(ky/2-1)/gamma(ky/2+1)/ky/n^2*PCvM
  
  out <- PCvM
  return(out)
}