flm_test <- function(B, n, X, pcX, pcY, hat_b)
{
  hat_y <- pcX$x%*%hat_b
  residuals <- pcY$x - hat_y
  
  ## WILD BOOTSTRAP RESAMPLING
  
  P = diag(1, nrow = n, ncol = n) - pcX$x %*% pseudoinverse(t(pcX$x) %*% pcX$x) %*% t(pcX$x)
  PCvM_star <- 0; 
  
  # Adot
  Adot.vec=Adot(X)
  
  # Obtain the entire matrix Adot
  Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
  Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
  Ad=t(Ad)
  Ad=Ad+t(Ad)-diag(diag(Ad))
  
  PCvM <- PCvM_statistic(fdata(pcX$x),residuals, Ad)
  
  for (i in 1:B) {
    #browser()
    res_star <- rwild(residuals, "golden")
    Y_star <- pcY$x - residuals + res_star
    res_hat_star = P %*% Y_star
    PCvM_star[i] = PCvM_statistic(fdata(pcX$x), res_hat_star, Ad)
  }
  
  hist(PCvM_star)
  pvalue = sum(PCvM_star > PCvM)/B
  return(pvalue)
}