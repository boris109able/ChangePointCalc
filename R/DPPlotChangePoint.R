DPPlotChangePoint <- function(sol, n, p, breaks=NULL)
{
  #sol is the solution by DP
  #n is the total number of data point
  #breaks is the break tick the user want to show on the plot, 
  #if breaks is NULL the default setting is to divide x-axis into 10 equal segment and label the 
  #tick accordingly
  t_sum = mat.or.vec(n, 1)
  for(i in 1:length(sol))
  {
    t_sum[sol[i]]<-1
  }
  t_sum[1]<-0
  t_sum[n]<-0
  plotData <- data.frame(pos=1:length(t_sum), normL1=t_sum)
  plotTsum <- ggplot(plotData, aes(pos, normL1))
  {
    if(is.null(breaks))
    {
      plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab("Change-point")+scale_x_continuous(breaks=seq(1, n, floor(n/10)))+scale_y_continuous(breaks=c(0,1))
    }
    else
    {
      plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab("Change-point")+scale_x_continuous(breaks=breaks)+scale_y_continuous(breaks=c(0,1))
    }
  }   
  #plotTsum + geom_point() + theme_bw(base_size=16) + xlab("Number of Observations") + ylab("Change-point")+scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000))+scale_y_continuous(breaks=c(0,1))
}