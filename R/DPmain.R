DPmain <- function(x, y, K)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  c1<-ComputeCost(x,y,n,p)
  sol <- DPSol(x,y,n,p,K,c1)#the fifth parametet is the number of intervals produced
  predictedAlpha <- PrintAlpha(x,y,sol)
  return(list(Sol=sol, AlphaHat=predictedAlpha))
}