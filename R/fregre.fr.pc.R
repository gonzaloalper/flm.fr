fregre.fr.pc=function (X, Y, k_x = 3, k_y = 3) 
  # Returns the regression parameters for a functional respone Y and a functional covariate X,
  # computed in two principal components bases, with k_y and k_x elements respectively.
{
  # Compute principal components
  pcX <- flm.pc(X, k_x) # k_x principal components of X
  pcY <- flm.pc(Y, k_y) # k_y principal components of Y
  
  # Dimensions
  x_col = ncol(pcX$x); y_col = ncol(pcY$x)
  
  # Linear regression using multivariate model and pseudo-inverse
  hat_beta <- t(pcX$x[,1:x_col])%*%pcX$x[,1:x_col]
  hat_beta <- pseudoinverse(hat_beta)
  hat_beta <- hat_beta%*%t(pcX$x[,1:x_col])
  hat_beta <- hat_beta%*%pcY$x[,1:y_col]

  out <- list(coefficients = hat_beta)
  class(out) = "fdata.comp"
  return(out)
}
