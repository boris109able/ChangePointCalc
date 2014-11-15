DPSol <- function(x,y,n,p,Kstar,c)
{
  #This function is the main function of DP
  #This function first computes the change-points via DP, then post-process the results into more
  #readable form
  result <- DPCalc(x,y,n,p,Kstar,c)
  r <- result$r
  s <- result$s
  K <- Kstar+1
  sol <- vector()
  sol <- append(sol,n)
  while(K>0)
  {
    sol <- append(sol,s[K,n])
    n <- s[K,n]
    K <- K-1
  }
  sol2 <- vector()
  len <- length(sol)
  for(i in len:1)
  {
    sol2 <- append(sol2,sol[i])
  }
  return(sol2)
}