RealData <- function(x, y)
{
  # data formatting using real data
  # x is a matrix n by p
  # y is a vector n by 1
  # each tuple(x[i], y[i]) is a data point
  n <- dim(x)[1]
  p <- dim(x)[2]
  XTilde <- matrix(0, n, np)
  for (i in 1:n)
  {
    for (j in 1:i)
    {
      XTilde[i, ((j-1)*p+1):(j*p)] <- x[i,]
    }
  }
  data <- list(x=XTilde, y=y)
  return(data)
}