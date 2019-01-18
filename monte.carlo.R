monte.carlo <- function(n, surface_beta, argvalsx, argvalsy, p, q, M, method = "FPCs")
{
  lx <- length(argvalsx); ly <- argvalsy;
  
  library(fda.usc)
  
  # Adot
  Adot.vec=Adot(X.fdata)
  
  # Obtain the entire matrix Adot
  Ad=diag(rep(Adot.vec[1],dim(X.fdata$data)[1]))
  Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
  Ad=t(Ad)
  Ad=Ad+t(Ad)-diag(diag(Ad))
  
  # Monte Carlo
  t <- list()
  pb <- txtProgressBar(style = 3)
  
  for (counter in 1:M) {
    # Data generation
    X.fdata <- rproc2fdata(n = n, argvalsx, sigma="brownian")
    # Error
    noise <- matrix(rnorm(n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)
    # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
    Y.fdata <- fdata(mdata = noise, argvals = argvalsy, rangeval = range(argvalsy), 
               names = list(main = "Functional response", ylab = "Y"))
    length_x <- diff(X$rangeval)
    
    for (j in 1:ly) {
      b <- surface_beta[, j]
      Y.fdata$data[, j] <- sapply(1:n, function(i) {
        integrateSimp1D(b * X.fdata$data[i, ], length_x)
      })
    }
    
    # Add noise
    d <- 0
    Y.fdata$data <- Y.fdata$data + noise + d * exp(X.fdata$data)
    
    n <- dim(X.fdata)[1]
    # Basis expansions
    X_range <- X.fdata$rangeval[2] - X$rangeval[1]
    Y_range <- Y.fdata$rangeval[2] - Y$rangeval[1]
    
    if (method = "FPCs")
    {
      X <- fpc(X.fdata,p,equispaced = TRUE)
      Y <- fpc(Y.fdata,q,equispaced = TRUE)
    }
    else if (method = "Fourier")
    {
      basis_fourier_x <- t(eval.basis(evalarg = X.fdata$argvals,
                                      basisobj = create.fourier.basis(rangeval = c(X.fdata$rangeval[1], X.fdata$rangeval[2]), 
                                                                      nbasis = p)))
      basis_fourier_y <- t(eval.basis(evalarg = Y.fdata$argvals,
                                      basisobj = create.fourier.basis(rangeval = c(Y.fdata$rangeval[1], Y.fdata$rangeval[2]), 
                                                                      nbasis = q)))
      for (i in 1:p){
        basis_fourier_x[i,] <-
          basis_fourier_x[i,]/sqrt(integrateSimp1D(basis_fourier_x[i,]^2, c(X_range)))
      }
      for (i in 1:q){
        basis_fourier_y[i,] <-
          basis_fourier_y[i,]/sqrt(integrateSimp1D(basis_fourier_y[i,]^2, c(Y_range)))
      }
      basis_fourier_x <- t(basis_fourier_x); basis_fourier_y <- t(basis_fourier_y)
      X <- matrix(0, nrow = n, ncol = p)
      Y <- matrix(0, nrow = n, ncol = q)
      
      for (i in 1:n){
        for (j in 1:p){
          projX[i,j] <- integrateSimp1D(X$data[i,]*basis_fourier_x[,j],X_range)
        }
      }
      for (i in 1:n){
        for (j in 1:q){
          projY[i,j] <- integrateSimp1D(Y$data[i,]*basis_fourier_y[,j],Y_range)
        }
      }
    }
    
    t[[counter]] <- flm.fr.test(projX, projY, B_boot <-  500, Ad)
    setTxtProgressBar(pb, counter / M)
    
  }
  
  # Density statistic vs densities bootstraped statistics
  par(mfrow = c(1,2))
  orig_stats <- sapply(t, function(x) x$orig_stat)
  plot(density(orig_stats), lwd = 5, xlim = c(0, max(orig_stats)))
  rug(orig_stats)
  sapply(t[1:50], function(x) lines(density(x$boot_stat), 
                                    col = rgb(1, 0, 0, alpha = 0.25)))
  
  # Distribution of p-values
  hist(sapply(t, function(x) x$p_value))
}