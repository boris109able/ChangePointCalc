DPCalc<-function(x,y,n,p,Kstar,c)
{
  #This function do the main calculation of computing change-points via DP
  r <- matrix(0,ncol=n,nrow=Kstar+1)
  s <- matrix(0,ncol=n,nrow=Kstar+1)
  for(t in 1:n)
  {
    r[1,t] <- c[1,t]
    s[1,t] <- 1
  }
  
  for(K in 2:(Kstar+1))
  {
    #print(K)
    for(t in K:n)
    {
      q=-1
      #print(t)
      for(j in (K-1):(t-1))
      {
        #print(j)
        #print(t)
        tmp <- r[K-1,j]+c[j+1,t]
        if(q>tmp || q<0)
        {
          q<-tmp
          s[K,t]<-j+1
        }       
      }
      r[K,t] <- q
    }
  }
  return(list(r=r,s=s))
}