SGLcalc <-function(data, n, p, thresh = 0.01, nlam = 20, lambdas = NULL, beta.naught = rep(0,ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 0.01, gamma = 0.8, step = 1, reset = 10, alpha = 0.95, min.frac = 0.05, verbose = FALSE)
{
  #thresh/n is the epsilon_innter
  #outer.thresh is the epsilon_outer
  X <- data$x
  y <- data$y
  n <- nrow(X)
  np <- ncol(X)
  #Prepocessing with C++ information
  #alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))
  nlam = length(lambdas)
  beta.old <- rep(0, p)
  beta.is.zero <- rep(1, n)
  beta <- array(0, c(ncol(X),nlam))
  eta <- rep(0,n)
  #Begin calculating for each lambda
  for(k in 1:nlam)
  {
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    eta <- rep(0,n)      
    #junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambdas[k]*alpha), lambda2 = as.double(lambdas[k]*(1-alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset),PACKAGE = "ChangePointCalc")
    junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), n = n, p = p, np = n*p, lambda1 = as.double(lambdas[k]*alpha), lambda2 = as.double(lambdas[k]*(1-alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset),PACKAGE = "ChangePointCalc")
   
    beta.new <- junk$beta
    beta[,k] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if(verbose == TRUE){
      write(paste("***Lambda", k, "***"),"")
    }
  }
  return(list(beta = beta[unOrd,], lambdas = lambdas))
}
