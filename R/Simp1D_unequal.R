Simp1D_unequal <- function (argvals, fx, na.rm = TRUE) 
{
  lengthInterval <- tail(argvals, n = 1) - argvals[1]
  n <- sum(!is.na(fx))
  h <- 0
  for (i in 1:n - 1){
    h[i] <- argvals[i+1] - argvals[i]
  }
  int <- h[1] * fx[1] + sum(h[2:n-2] * fx[2:n-2], na.rm = na.rm) + h[n - 1] * fx[n - 1]
  int <- h * int
  return(int)
}