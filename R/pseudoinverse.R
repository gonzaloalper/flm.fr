pseudoinverse = function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  inv <- try(solve(X), silent = TRUE) # try to compute the inverse
  
  if (class(inv) != "try-error") {
    inv
  }
  else {
    X_svd <- svd(X) # singular value decomposition
    if (is.complex(X)) 
      X_svd$u <- Conj(X_svd$u) #  left singular vectors of X
    d_p <- X_svd$d > max(tol * X_svd$d[1L], 0) # posotive singular values of X
    warning("System is computationally singular (rank  ", 
            dim(X)[2], ")\n\n          The  matrix inverse is computed by svd (effective rank ", 
            sum(d_p), ")")
    if (all(d_p)) 
      X_svd$v %*% (1/X_svd$d * t(X_svd$u)) 
    else if (!any(d_p)) 
      array(0, dim(X)[2L:1L])
    else X_svd$v[, d_p, drop = FALSE] %*% ((1/X_svd$d[d_p]) * 
                                                 t(X_svd$u[, d_p, drop = FALSE]))
  }
}