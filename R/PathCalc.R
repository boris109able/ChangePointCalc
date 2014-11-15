PathCalc <- function(data, index, alpha = 0.95, min.frac = 0.05, nlam = 20){

reset <- 10
step <- 1
gamma <- 0.8

inner.iter <- 1000
outer.iter <- 1000
thresh = 10^(-2)
outer.thresh = thresh
  
  n <- nrow(data$x)

  X <- data$x
  resp <- data$y
  n <- nrow(X)
  p <- ncol(X)

  ## Setting up group lasso stuff ##
     
  
  ## Coming up with other C++ info ##
    
  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0,(num.groups+1))
  for(i in 1:num.groups){
    range.group.ind[i] <- min(which(index == groups[i])) - 1
  }
  range.group.ind[num.groups + 1] <- ncol(X)
  
  group.length <- diff(range.group.ind)
  

lambda.max <- rep(0,num.groups)

if((alpha != 0)*(alpha != 1)){

for(i in 1:num.groups){
  ind <- groups[i]
  X.fit <- X[,which(index == ind)]
  cors <- t(X.fit) %*% resp
  ord.cors <- sort(abs(cors), decreasing = TRUE)

  if(length(ord.cors) > 1){ 
  norms <- rep(0,length(cors)-1)

  lam <- ord.cors/alpha


  for(j in 1:(length(ord.cors)-1)){
    norms[j] <- sqrt(sum((ord.cors[1:j]-ord.cors[j+1])^2))
  }
  if(norms[1] >= lam[2] * (1-alpha)*sqrt(group.length[i])){
    our.cors <- ord.cors[1]
    our.range <- c(ord.cors[2], ord.cors[1])/alpha
  }else{
    if(norms[length(ord.cors)-1] <= lam[length(ord.cors)] * (1-alpha)*sqrt(group.length[i])){
      our.cors <- ord.cors
      our.range <- c(0, ord.cors[length(ord.cors)])/alpha
    } else{
      my.ind <- max(which(norms[-length(norms)] <= lam[2:(length(norms))] * (1-alpha) * sqrt(group.length[i]))) + 1
      our.cors <- ord.cors[1:my.ind]
      our.range <- c(ord.cors[my.ind+1], ord.cors[my.ind])/alpha
    }
  }
  nn <- length(our.cors)
  if(alpha == 0.5){
    alpha = 0.500001
  }
  
  A.term <- nn*alpha^2 - (1 - alpha)^2*group.length[i]
  B.term <- - 2 * alpha * sum(our.cors)
  C.term <- sum(our.cors^2)

  lams <- c((-B.term + sqrt(B.term^2 - 4 * A.term * C.term))/(2*A.term), (-B.term - sqrt(B.term^2 - 4 * A.term * C.term))/(2*A.term))
  
  
  lambda.max[i] <- min(subset(lams, lams >= our.range[1] & lams <= our.range[2]))
  }
     if(length(ord.cors) == 1){
       lambda.max[i] <- ord.cors
     }
}
}
if(alpha == 1){
  lambda.max <- abs(t(X) %*% resp)
}
if(alpha == 0){
  for(i in 1:num.groups){
  ind <- groups[i]
  X.fit <- X[,which(index == ind)]
  cors <- t(X.fit) %*% resp
  lambda.max[i] <- sqrt(sum(cors^2)) / sqrt(group.length[i])
  }
}
max.lam <- max(lambda.max)
min.lam <- min.frac*max.lam
lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))
return(lambdas/nrow(X))
}
