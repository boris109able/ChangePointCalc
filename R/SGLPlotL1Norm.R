SGLPlotL1Norm <- function(fit, n, p, num=1, breaks=NULL)
{
  # this function plot the l1 norm of theta in the result
  # num is the num th lambda we used
  # n is the number of data sample
  # p is the number of dimension
  #breaks is the break tick the user want to show on the plot, 
  #if breaks is NULL the default setting is to divide x-axis into 10 equal segment and label the 
  #tick accordingly
  t_sum = matrix(NA, n, 1)
  for (i in 1:n)
  { t_sum[i] = sum(abs(fit$beta[((i-1)*p+1):(i*p), num]))}
  plotData <- data.frame(pos=1:length(t_sum), normL1=t_sum)
  plotTsum <- ggplot(plotData, aes(pos, normL1))
  {
    if(is.null(breaks))
    {
      plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab(expression(paste(group("|",group("|", theta[i], "|"),"|"))[1]))+scale_x_continuous(breaks=seq(1, n, floor(n/10)))
    }
    else
    {
      plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab(expression(paste(group("|",group("|", theta[i], "|"),"|"))[1]))+scale_x_continuous(breaks=breaks)
    }
  }
  #plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab(expression(paste(group("|",group("|", theta[i], "|"),"|"))[1]))+scale_x_continuous(breaks=seq(1, n, floor(n/10)))
}
