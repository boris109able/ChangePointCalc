SGLLogitmain <- function(data, gamma = 0.9, maxit = 2000, thresh = 0.0005, min.frac = 0.05, nlam = 20, alpha = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, lambdas = NULL)
{
  n <- nrow(data$x)
  p <- ncol(data$x)
  np <- n*p
  
  X.transform <- NULL
  #standardize = FALSE
  if(standardize == TRUE){#this part has error since some of var can be zero!!! revise it below
    X <- data$x
    #means <- apply(X,2,mean)
    #X <- t(t(X) - means)
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
    X.transform <- var
  }
  
  #if(standardize == TRUE){
  #  intercept <- mean(data$y)
  #  data$y <- data$y - intercept
  #}  
  
  # if lambdas is not provided, then calculate lambdas
  if(is.null(lambdas)){
    lambdas <- PathCalcLogistic(data = data, gamma=gamma, min.frac = min.frac, nlam = nlam)
  }
  
  
  Sol <- SGLLogitcalc(data, thresh=thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)
  if(standardize == TRUE){
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas,  X.transform = X.transform)
  }
  if(standardize == FALSE){
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas,  X.transform = NULL)
  }
  class(Sol) = "SGLLogit"
  return(Sol)
}


