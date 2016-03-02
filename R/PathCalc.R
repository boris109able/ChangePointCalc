PathCalc <- function(data, gamma = 0.90, min.frac = 0.05, nlam = 20){
 
  X <- data$x
  resp <- data$y
  n <- nrow(X)
  p <- ncol(X)
  np <- n*p
  num.groups <- n

lambda.max <- rep(0,num.groups)

if((gamma != 0)*(gamma != 1)){

for(i in 1:num.groups){
  if(i<num.groups)
  {
    cors <- t(X[i:n,]) %*% resp[i:n]*2/n
  }
  else {
    cors <- t(X[i:n,] %*% t(resp[i:n]))*2/n
  }
  ord.cors <- sort(abs(cors), decreasing = TRUE)

  if(length(ord.cors) > 1){ 
  norms <- rep(0,length(cors)-1)
  lam <- ord.cors/(1-gamma)


  for(j in 1:(length(ord.cors)-1)){
    norms[j] <- sqrt(sum((ord.cors[1:j]-ord.cors[j+1])^2))
  }
  if(norms[1] >= lam[2] * gamma){
    our.cors <- ord.cors[1]
    our.range <- c(ord.cors[2], ord.cors[1])/(1-gamma)
  }else{
    if(norms[length(ord.cors)-1] <= lam[length(ord.cors)] * gamma){
      our.cors <- ord.cors
      our.range <- c(0, ord.cors[length(ord.cors)])/(1-gamma)
    } else{
      my.ind <- max(which(norms[-length(norms)] <= lam[2:(length(norms))] * gamma)) + 1
      our.cors <- ord.cors[1:my.ind]
      our.range <- c(ord.cors[my.ind+1], ord.cors[my.ind])/(1-gamma)
    }
  }
  nn <- length(our.cors)
  if(gamma == 0.5){
    gamma = 0.500001
  }
  
  A.term <- nn*(1-gamma)^2 - gamma^2
  B.term <- - 2 * (1-gamma) * sum(our.cors)
  C.term <- sum(our.cors^2)
  lams <- c((-B.term + sqrt(B.term^2 - 4 * A.term * C.term))/(2*A.term), (-B.term - sqrt(B.term^2 - 4 * A.term * C.term))/(2*A.term))
  if((lams[1] >= our.range[1] && lams[1] <= our.range[2]) || (lams[2] >= our.range[1] && lams[2] <= our.range[2])) {
    lambda.max[i] <- min(subset(lams, lams >= our.range[1] & lams <= our.range[2]))
  }
  else {
    lambda.max[i] = our.range[1]
  }
  }
     if(length(ord.cors) == 1){
       lambda.max[i] <- ord.cors
     }
}
}
if(gamma == 0){
  lambda.max <- abs(t(X) %*% resp)
}
if(gamma == 1){
  for(i in 1:num.groups){
    if (i<num.groups) {
    cors <- t(X[i:n,]) %*% resp[i:n] * 2/n}
    else {
      cors <- t(X[i:n,] %*% t(resp[i:n]))*2/n
    }
  lambda.max[i] <- sqrt(sum(cors^2))
  }
}
max.lam <- max(lambda.max)
min.lam <- min.frac*max.lam
lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))
return(lambdas)
}
