integrateSimp2D <- function(fxy, lengthInterval = rep(2 * pi, 2),
                            na.rm = TRUE) {
  
  # Length and spacing of the grid
  n <- dim(fxy)
  h <- lengthInterval / (n - 1)
  
  # Iteration of equation (4.1.14) in Numerical Recipes in F77
  # Integral on X
  intX <- 3 / 8 * (fxy[1, ] + fxy[n[1], ]) +
    7 / 6 * (fxy[2, ] + fxy[n[1] - 1, ]) +
    23 / 24 * (fxy[3, ] + fxy[n[1] - 2, ]) +
    colSums(fxy[4:(n[1] - 3), ], na.rm = na.rm)
  
  # Integral on Y
  int <- 3 / 8 * (intX[1] + intX[n[2]]) +
    7 / 6 * (intX[2] + intX[n[2] - 1]) +
    23 / 24 * (intX[3] + intX[n[2] - 2]) +
    sum(intX[4:(n[2] - 3)], na.rm = na.rm)
  int <- int * prod(h)
  
  return(int)
  
}
