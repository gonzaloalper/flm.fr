monte.carlo <- function(n, surface_beta, argvalsx, argvalsy, p, q, M, method = "FPCs")
{
  
  lx <- length(argvalsx)
  ly <- length(argvalsy) # TODO: ly estaba como ly <- argvalsy
  
  # Adot
  Adot.vec <- fda.usc::Adot(X.fdata)
  
  # Obtain the entire matrix Adot
  Ad <- diag(rep(Adot.vec[1], dim(X.fdata$data)[1]))
  Ad[upper.tri(Ad, diag = FALSE)] <- Adot.vec[-1]
  Ad <- t(Ad)
  Ad <- Ad + t(Ad) - diag(diag(Ad))
  
  # Monte Carlo
  t <- list()
  pb <- txtProgressBar(style = 3)
  
  for (counter in 1:M) {
    
    # Data generation
    X.fdata <- fda.usc::rproc2fdata(n = n, argvalsx, sigma="brownian")
    # Error
    noise <- matrix(rnorm(n * ly, mean = 0, sd = 0.1), nrow = n, ncol = ly)
    # Response: Y(t) = \int X(s) beta(s, t) ds + eps(t)
    Y.fdata <- fda.usc::fdata(mdata = noise, argvals = argvalsy, 
                              rangeval = range(argvalsy), 
                              names = list(main = "Functional response", 
                                           ylab = "Y"))
    length_x <- diff(X$rangeval)
    
    for (j in 1:ly) {
      b <- surface_beta[, j]
      Y.fdata$data[, j] <- sapply(1:n, function(i) {
        sdetorus::integrateSimp1D(b * X.fdata$data[i, ], length_x)
      })
    }
    
    # Add noise
    d <- 0
    Y.fdata$data <- Y.fdata$data + noise + d * exp(X.fdata$data)
    
    n <- dim(X.fdata)[1]
    # Basis expansions
    X_range <- diff(X.fdata$rangeval)
    Y_range <- diff(Y.fdata$rangeval)
    
    if (method = "FPCs") {
      
      X <- fpc(X.fdata, p, equispaced = TRUE)
      Y <- fpc(Y.fdata, q, equispaced = TRUE)
      
    } else if (method = "Fourier") {
      
      basis_fourier_x <- t(fda::eval.basis(evalarg = X.fdata$argvals,
                                           basisobj = 
                                             fda::create.fourier.basis(
                                               rangeval = X.fdata$rangeval, 
                                               nbasis = p)))
      basis_fourier_y <- t(fda::eval.basis(evalarg = Y.fdata$argvals,
                                           basisobj = 
                                             fda::create.fourier.basis(
                                               rangeval = Y.fdata$rangeval, 
                                               nbasis = q)))
      for (i in 1:p) {
        basis_fourier_x[i, ] <- basis_fourier_x[i, ] / 
          sqrt(sdetorus::integrateSimp1D(basis_fourier_x[i, ]^2, X_range))
      }
      for (i in 1:q) {
        basis_fourier_y[i, ] <- basis_fourier_y[i, ] /
          sqrt(sdetorus::integrateSimp1D(basis_fourier_y[i, ]^2, Y_range))
      }
      basis_fourier_x <- t(basis_fourier_x) # TODO: no tiene sentido hacer dos veces traspuesta
      basis_fourier_y <- t(basis_fourier_y)
      X <- matrix(0, nrow = n, ncol = p)
      Y <- matrix(0, nrow = n, ncol = q)
      
      for (i in 1:n){
        for (j in 1:p){
          projX[i, j] <- sdetorus::integrateSimp1D(X$data[i, ] * basis_fourier_x[, j], X_range)
        }
      }
      for (i in 1:n){
        for (j in 1:q){
          projY[i, j] <- sdetorus::integrateSimp1D(Y$data[i,] * basis_fourier_y[, j], Y_range)
        }
      }
    }
    
    t[[counter]] <- flm.fr.test(projX, projY, B_boot <-  500, Ad)
    setTxtProgressBar(pb, counter / M)
    
  }
  
  # Density statistic vs densities bootstraped statistics
  par(mfrow = c(1, 2))
  orig_stats <- sapply(t, function(x) x$orig_stat)
  plot(density(orig_stats), lwd = 5, xlim = c(0, max(orig_stats)))
  rug(orig_stats)
  sapply(t[1:50], function(x) lines(density(x$boot_stat), 
                                    col = rgb(1, 0, 0, alpha = 0.25)))
  
  # Distribution of p-values
  hist(sapply(t, function(x) x$p_value))
  
}