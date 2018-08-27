PCvM_statistic = function(X, residuals, Ad)
{
  n <- dim(residuals)[1]
  ky <- dim(residuals)[2]
  
  PCvM <- sum(diag(t(residuals) %*% Ad %*% residuals)) * pi^(ky/2-1)/gamma(ky/2+1)/ky/n^2
  
  out <- PCvM
  return(out)
}
