fourier_expansion = function (X, Y, n_basis_x, n_basis_y, lx, ly) 
{
  X_range <- X$rangeval[2] - X$rangeval[1]
  Y_range <- Y$rangeval[2] - Y$rangeval[1]
  
  basis_fourier_x <- t(eval.basis(evalarg = seq(X$rangeval[1], X$rangeval[2], l = lx),
                                  basisobj = create.fourier.basis(rangeval = c(X$rangeval[1], X$rangeval[2]), 
                                                                  nbasis = n_basis_x)))
  basis_fourier_y <- t(eval.basis(evalarg = seq(Y$rangeval[1], Y$rangeval[2], l = ly),
                                  basisobj = create.fourier.basis(rangeval = c(Y$rangeval[1], Y$rangeval[2]), 
                                                                  nbasis = n_basis_y)))
  
  # Check that norm of the basis for X
  #for (i in 1:n_basis_x){
  #  print(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
  #}
  
  # Normalize the basis functions for X
  for (i in 1:n_basis_x){
    basis_fourier_x[i,] <-
      basis_fourier_x[i,]/sqrt(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
  }
  
  # Check that the basis or X has functional norm 1
  #for (i in 1:n_basis_x){
  #  print(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
  #}
  
  # Normalize the basis functions for Y
  for (i in 1:n_basis_y){
    basis_fourier_y[i,] <-
      basis_fourier_y[i,]/sqrt(integrateSimp1D(basis_fourier_y[i,]^2, c(Y_range)))
  }
  
  basis_fourier_x <- t(basis_fourier_x); basis_fourier_y <- t(basis_fourier_y)
  
  # par(mfrow = c(1, 1)); matplot(basis_fourier_x, type="l")
  
  projX <- matrix(0, nrow = n, ncol = n_basis_x)
  projY <- matrix(0, nrow = n, ncol = n_basis_y)
  
  # Fiind the projections
  for (i in 1:n){
    for (j in 1:n_basis_x){
      projX[i,j] <- integrateSimp1D(X$data[i,]*basis_fourier_x[,j],X_range)
    }
  }
  
  for (i in 1:n){
    for (j in 1:n_basis_y){
      projY[i,j] <- integrateSimp1D(Y$data[i,]*basis_fourier_y[,j],Y_range)
    }
  }
  
  out <- list(projX, projY, basis_fourier_x, basis_fourier_y)
  return(out)
}
