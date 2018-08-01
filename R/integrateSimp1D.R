integrateSimp1D <- function(fx, lengthInterval = 2 * pi, na.rm = TRUE) {
  
  # Length and spacing of the grid
  n <- sum(!is.na(fx))
  h <- lengthInterval / (n - 1)
  
  # Equation (4.1.14) in Numerical Recipes in F77.
  # Much better than equation (4.1.13).
  int <- 3 / 8 * (fx[1] + fx[n]) +
    7 / 6 * (fx[2] + fx[n - 1]) +
    23 / 24 * (fx[3] + fx[n - 2]) +
    sum(fx[4:(n - 3)], na.rm = na.rm)
  int <- h * int
  
  return(int)
  
}
