beta_hat_net <- function(x, y, type = c("ridge", "lasso")[1], lambda = NULL, 
                         k_folds = nrow(x), ...) {
  
  if (is.null(dim(x)) | is.null(dim(x)) | ncol(x) == 1 | ncol(y) == 1) {
    
    stop("x and y must be matrices of two or more columns")
    
  }
  
  alpha <- switch(type, "ridge" = 0, "lasso" = 1)
  if (is.null(lambda)) {
    
    # Ensure that the optimal lambda does not fall in one extreme
    keep <- TRUE
    while (keep) {
      
      cv <- glmnet::cv.glmnet(x = x, y = y, family = "mgaussian", alpha = alpha, 
                              lambda = lambda, nfolds = k_folds, 
                              intercept = FALSE, standardize = FALSE, 
                              standardize.response = FALSE, ...)
      
      # Expand grid if optimum is in the 10% lower part of the grid or in the 
      # 90% upper part of the grid
      
      if (abs(which.min(cv$cvm) / length(cv$lambda) - 0.5) > 0.45) {
        
        if (cv$lambda[which.min(cv$cvm)] < 1e-4) {
          
          keep <- FALSE
          
        } else {
          
          message(paste0("lambda.min = ", sprintf("%.6f", cv$lambda.min), 
                         " close to boundary (", sprintf("%.6f", min(cv$lambda)), 
                         ", ", sprintf("%.6f", max(cv$lambda)), 
                         "), expanding search grid"))
          lower_lambda <- cv$lambda[which.min(cv$cvm)]/4
          upper_lambda <- cv$lambda[which.min(cv$cvm)] + 
            cv$lambda[which.min(cv$cvm)]/2
          lambda <- seq(lower_lambda, upper_lambda, 
                        length.out = length(cv$lambda))
          lambda[which(lambda<0)] <- 10^lambda[which(lambda<0)]
          
        }
        
      } else {
        
        message(paste0("lambda.min = ", sprintf("%.6f", cv$lambda.min), 
                       " converging to zero, terminating search"))
        keep <- FALSE
        
      }
      
    }
    
    lambda <- cv$lambda.1se
    
    beta_hat <- predict(cv, type = "coefficients", s = lambda)
    beta_hat <- do.call(cbind, args = as.list(beta_hat))[-1, ]
    
  } else {
    
    fit <- glmnet::glmnet(x = x, y = y, family = "mgaussian", alpha = alpha, 
                          lambda = lambda, intercept = FALSE, 
                          standardize = TRUE, standardize.response = FALSE, ...)
    beta_hat <- do.call(cbind, args = as.list(fit$beta))
    
  }
  
  return(list("beta_hat" = beta_hat, "lambda" = lambda))
  
}