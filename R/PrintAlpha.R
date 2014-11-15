PrintAlpha<-function(x,y,sol)
{
  #This function return the alpha solution
  n<-dim(x)[1]
  p<-dim(x)[2]
  K<-length(sol)
  Alpha <- matrix(NA,K-1,p)
  if(K==2)
  {
    Alpha[i,]<-CalcAlpha(x,y,n,p,sol[1],sol[2])
  }
  else
  {
    for(i in 1:(K-1))
    {
      if(i==(K-1))
      {
        Alpha[i,]<-CalcAlpha(x,y,n,p,sol[i],sol[i+1])
      }
      else
      {
        Alpha[i,]<-CalcAlpha(x,y,n,p,sol[i],sol[i+1]-1)
      }
    }
  }
  return(Alpha)
}