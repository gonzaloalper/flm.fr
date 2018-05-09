reconstruct = function (fdataobj, pcX) 
{
  par(mfrow = c(1, 2))
  
  # Plot original data
  plot(fdataobj[1,])
  lines(fdataobj[2,],col=2)
  lines(fdataobj[3,],col=3)
  
  # Plot reconstructed data with PC's
  XReconstructed <- pcX$x%*%t(pcX$rotation)
  plot(XReconstructed[1,],type="l")
  lines(XReconstructed[2,],col=2)
  lines(XReconstructed[3,],col=3)
}
