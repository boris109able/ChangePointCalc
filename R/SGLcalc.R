SGLcalc <-function(data, gamma = 0.9, thresh = 0.01, nlam = 20, lambdas = NULL, beta.naught = rep(0,ncol(data$x)), inner.iter = 1000, outer.iter = 1000, outer.thresh = 0.01, alpha = 0.8, step = 1, reset = 10,  min.frac = 0.05, verbose = FALSE)
{
  #thresh/n is the epsilon_innter
  #outer.thresh is the epsilon_outer
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)
  np <- n*p
  #Prepocessing with C++ information
  #alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))
  nlam = length(lambdas)
  beta.old <- rep(0, np)
  beta.is.zero <- rep(1, n)
  beta <- array(0, c(np,nlam))
  eta <- rep(0,n)
  #Begin calculating for each lambda
  for(k in 1:nlam)
  {
    beta.is.zero <- rep(1, n)
    beta.old <- rep(0, np)
    eta <- rep(0,n)      
    #junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambdas[k]*alpha), lambda2 = as.double(lambdas[k]*(1-alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset),PACKAGE = "ChangePointCalc")
    junk <- .C("solveSGL", X = as.double(as.vector(X)), y = as.double(y), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), lambda1 = as.double(lambdas[k]*(1-gamma)), lambda2 = as.double(lambdas[k]*gamma), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(alpha), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset),PACKAGE = "ChangePointCalc")
   
    beta.new <- junk$beta
    beta[,k] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if(verbose == TRUE){
      write(paste("***Lambda", k, "***"),"")
    }
  }
  return(list(beta = beta, lambdas = lambdas))
}
