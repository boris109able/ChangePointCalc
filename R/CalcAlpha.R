CalcAlpha<-function(x,y,n,p,i,j)
{
  A<-matrix(0,p,p)
  B<-matrix(0,p,1)
  for(t in i:j)
  {
    A <- A + x[t,]%*%t(x[t,])
    B <- B + y[t]*x[t,]
  }
  alpha <- ginv(A)%*%B
  return(alpha) 
}