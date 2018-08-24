fregre.pc.ex <- function(B)
{
  library(fda.usc)
  library(viridis)
  # Execute all the functions in the R directory from the above file:
  #sourceDir <- function(directory = "R", ...) {
  #  invisible(sapply(dir(directory),
  #                   function(x) source(paste0(directory, "/", x), ...)))
  #}
  #sourceDir()
  
  # Argvals
  lx <- 205
  ly <- 205
  ss <- seq(0, 1, l = lx) # Argvals of X
  tt <- seq(0, 1, l = ly) # Argvals of Y
  
  # Sample covariate
  n <- 100
  X <- rproc2fdata(n = n, ss, sigma="brownian")
  #X <- r.ou(n = n, x0 = seq(-10, 10, l = n), alpha = 2, t = ss)
  X_range <- X$rangeval[2] - X$rangeval[1]
  hx <- (X$rangeval[2]-X$rangeval[1])/length(ss)
  
  # Theoretical beta
  beta <- function(s, t) {
    #diag(2,lx,ly)
    #2+(s - s + t - t)
    #cos(2 * pi * t * s) / (1 + (s - 0.5)^2) 
    pi^2*(s^2-1/(0.25+t)) * tanh(1+s + t^3)
    #cos(t * s / pi)^2
  }
  
  # Visualization
  surface_beta <- outer(ss, tt, FUN = beta)
  
  Y <- linear_model(X, tt, 0.1, beta)
  
  Y_range <- Y$rangeval[2] - Y$rangeval[1]
  
  # Plots
  par(mfrow = c(1, 2))
  plot(X)
  plot(Y)
  
  npcX<-25
  npcY<-25
  pcX <- fpc(X,npcX,equispaced = TRUE)
  pcY <- fpc(Y,npcY,equispaced = TRUE)
  
  hat_b <- t(pcX$x[,1:npcX])%*%pcX$x[,1:npcX]
  hat_b <- pseudoinverse(hat_b)
  hat_b <- hat_b%*%t(pcX$x[,1:npcX])
  hat_b <- hat_b%*%pcY$x[,1:npcY]
  
  theoretical_beta <- matrix(0, nrow = npcX, ncol = npcY)
  
  for (i in 1:npcX){
    for (j in 1:npcY){
      outer_product <- outer(pcX$rotation[,i],pcY$rotation[,j])
      acc <- surface_beta * outer_product
      theoretical_beta[i,j] <- integrateSimp2D(acc, c(X_range, Y_range))
    }
  }
  
  #EXPLAINED VARIANCE
  
  diferencia <- (sum(((hat_b - theoretical_beta)/theoretical_beta)^2))^(1/(npcX*npcY)); 
  explained_variance_X <- sum(pcX$d[1:npcX])/sum(pcX$d); explained_variance_X
  explained_variance_Y <- sum(pcY$d[1:npcY])/sum(pcY$d); explained_variance_Y
  
  ## SURFACES
  surface_beta_hat <- matrix(0, nrow = lx, ncol = ly)
  surface_pca <- matrix(0, nrow = lx, ncol = ly)
  for (i in 1:npcX){
    for (j in 1:npcY){
      outer_product <- outer(pcX$rotation[,i], pcY$rotation[,j])
      acc1 <- theoretical_beta[i,j] * outer_product
      surface_pca <- surface_pca + acc1
      acc2 <- hat_b[i,j] * outer_product
      surface_beta_hat <- surface_beta_hat + acc2
    }
  }
  
  par(mfrow = c(1, 3))
  minima <- c(min(surface_beta), min(surface_pca), min(surface_beta_hat))
  maxima <- c(max(surface_beta), max(surface_pca), max(surface_beta_hat))
  limits = c(min(minima), max(maxima))
  image(ss, tt, surface_beta, col = viridis(20), zlim = limits)
  image(ss, tt, surface_pca, col = viridis(20), zlim = limits)
  image(ss, tt, surface_beta_hat, col = viridis(20), zlim = limits)
  
  ## RESIDUALS AND BOOTSTRAP
  y_hat <- Y
  y_hat$data <- pcY$x %*% t(pcY$rotation)
  fresiduals <- Y - y_hat$data
  res_star <- fresiduals; for (j in 1:dim(fresiduals$data)[1]){ #PERTURBED RESIDUALS
    res_star$data[j,] <- fresiduals$data[j,] * rwild(1, "golden")}
  Y_star <- y_hat - fresiduals + res_star
  
  par(mfrow = c(1, 3))
  mm <- c(min(as.vector(fresiduals$data)), min(as.vector(res_star$data)))
  mx <- c(max(as.vector(fresiduals$data)), max(as.vector(res_star$data)))
  lims <- c(min(mm), max(mx))
  boxplot(as.vector(fresiduals$data), ylim = lims)
  boxplot(as.vector(res_star$data), ylim = lims)
  plot(ecdf(as.vector(fresiduals$data)),col=4)
  lines(ecdf(as.vector(res_star$data)),col=2)
  legend('topleft', legend=c("Original residuals", "Perturbed residuals"),col=c(3,2), lty=1, cex=0.8)
  
  summary(as.vector(fresiduals$data))
  summary(as.vector(res_star$data))
  
  # Proportion of the residuals that are larger than the bootstraped
  100*sum(fresiduals$data>res_star$data)/dim(fresiduals$data)[1]/dim(fresiduals$data)[2]
  
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
  
  residuals <- P %*% Y$data %*% pcY$rotation
  PCvM <- PCvM_statistic(pcX$x, residuals, Ad)
  
  Y_star <- Y
  res_hat_star <- res_star
  
  for (i in 1:B) {
    for (j in 1:dim(fresiduals$data)[1]){
      res_star$data[j,] <- fresiduals$data[j,] * rwild(1, "golden")
    }
    Y_star <- y_hat - fresiduals + res_star
    pcY_star <- Y_star$data %*% pcY$rotation
    res_hat_star <- P %*% pcY_star
    PCvM_star[i] = PCvM_statistic(pcX$x, res_hat_star, Ad)
  }
  
  par(mfrow = c(1, 1))
  hist(PCvM_star)
  abline(v = PCvM, col = 2)
  summary(PCvM_star)
  #print(PCvM)
  
  pvalue = sum(PCvM_star > PCvM)/B; print(pvalue)
  
  par(mfrow = c(2, 3))
  plot(Y)
  plot(y_hat)
  plot(Y_star)
  plot(fresiduals)
  plot(res_star)
  
  return(pvalue)
}
