fregre.lr.pc=function (fdataobj, Y, ncomp = 3) 
{
  # Compute principal components
  pcX <- flm.pc(fdataobj)
  pcY <- flm.pc(Y)
  
  # Dimensions
  x_row = nrow(pcX$x); x_col = ncol(pcX$x)
  y_row = nrow(pcY$x); x_col = ncol(pcY$x)
  
  # Linear regression using multivariate model and pseudo-inverse
  hat_beta <- t(pcX$x[,1:x_col])%*%pcX$x[,1:x_col]
  hat_beta <- pseudoinverse(hat_beta)
  hat_beta <- hat_beta%*%t(pcX$x[,1:x_col])
  hat_beta <- hat_beta%*%pcY$x[,1:y_col]

  out <- list(coefficients = hat_beta)
  class(out) = "fdata.comp"
  return(out)
}
