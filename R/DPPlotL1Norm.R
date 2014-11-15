DPPlotL1Norm<-function(predictedAlpha, sol, n, p, breaks=NULL)
{
  #predictedAlpha is the alpha calculated by DP
  #sol is the solution by DP
  #n is the total number of data point
  #p is the dimension of each data point
  #breaks is the break tick the user want to show on the plot, 
  #if breaks is NULL the default setting is to divide x-axis into 10 equal segment and label the 
  #tick accordingly
  #p <- dim(predictedAlpha)[2]
  predictedTheta <- matrix(0,length(sol)-1,p)#calculates theta
  for(i in 1:(length(sol)-1))
  {
    if(i==1)
    {
      predictedTheta[i,] <- predictedAlpha[i,]
    }
    else
    {
      predictedTheta[i,] <- predictedAlpha[i,] - predictedAlpha[i-1,]
    } 
  }
  #plotting
  t_sum = mat.or.vec(n, 1)
  for(i in 1:(length(sol)-1))
  {
    t_sum[sol[i]]<-sum(abs(predictedTheta[i,]))
  }
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
  #plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab(expression(paste(group("|",group("|", theta[i], "|"),"|"))[1]))+scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000))
}