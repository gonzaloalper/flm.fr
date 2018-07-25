trap1D_unequal <- function (argvals, fx, na.rm = TRUE) 
{
  lengthInterval <- tail(argvals, n = 1) - argvals[1]
  n <- sum(!is.na(fx))
  h <- diff(argvals)
  int <-   h[2:length(h)] * 0.5 * (fx[1:length(h)-1] - fx[2:length(h)])
  return(int)
}
