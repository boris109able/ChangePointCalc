 SGLmain <- function(data, maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL)
{
   n <- nrow(data$x)
   np <- ncol(data$x)
   p <- np/n

  X.transform <- NULL

  if(standardize == TRUE){#this part has error since some of var can be zero!!! revise it below
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    for (i in 1:length(var))
    {
      if (var[i]==0)#revised
      {
        var[i] <- 1
      }
    }
    X <- t(t(X) / var)
    data$x <- X
    X.transform <- list(X.scale = var, X.means = means)
  }
  # if lambdas is not provided, then calculate lambdas
  if(is.null(lambdas)){
    lambdas <- PathCalc(data = data, index=ceiling(1:np/p), alpha=alpha, min.frac = min.frac, nlam = nlam)
  }
  
  if(standardize == TRUE){
    intercept <- mean(data$y)
    data$y <- data$y - intercept
  }   
  
  Sol <- SGLcalc(data, thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)
  if(standardize == TRUE){
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, intercept = intercept, X.transform = X.transform)
  }
  if(standardize == FALSE){
     Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, intercept = NULL, X.transform = NULL)
  }
  class(Sol) = "SGL"
  return(Sol)
}


