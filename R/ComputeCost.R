ComputeCost<-function(x,y,n,p)
{
  #This function computes the cost of any interval [i,j]
  #The input is the data
  #The output is the cost matrix
  c <- matrix(0,n,n)
  for(i in 1:n)
  {
    #print(i)
    A <- matrix(0,p,p)
    B <- matrix(0,p,1)
    for(j in i:n)
    {
      #cat(sprintf("\"%d\" \"%d\"\n", i, j))
      A <- A + x[j,]%*%t(x[j,])
      B <- B + y[j]*x[j,] 
      alpha <- ginv(A)%*%B
      tmp <- 0
      for(t in i:j)
      {
        tmp <- tmp + (y[t]-(t(x[t,])%*%alpha)[1])^2
      }
      c[i,j] <- tmp
    }
  }
  return(c)
}