\name{DPPlotChangePoint}
\alias{DPPlotChangePoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot change-points indicated by 1 or 0 using DP approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot change-points indicated by 1 or 0 using DP approach
}
\usage{
DPPlotChangePoint(sol, n, p, breaks = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sol}{
%%     ~~Describe \code{sol} here~~
solution of DP approach
}
  \item{n}{
%%     ~~Describe \code{n} here~~
total number of data
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
list of breaks where to put a tick on the x-axis
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Bingwen Zhang, Jun Geng and Lifeng Lai, Detecting changes in regression models via sparse group lasso,
http://users.wpi.edu/~bzhang/icassp_v1.pdf
}
\author{
%%  ~~who you are~~
Bingwen Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sol, n, p, breaks = NULL) 
{
    t_sum = mat.or.vec(n, 1)
    for (i in 1:length(sol)) {
        t_sum[sol[i]] <- 1
    }
    t_sum[1] <- 0
    t_sum[n] <- 0
    plotData <- data.frame(pos = 1:length(t_sum), normL1 = t_sum)
    plotTsum <- ggplot(plotData, aes(pos, normL1))
    {
        if (is.null(breaks)) {
            plotTsum + geom_point() + theme_bw(base_size = 16) + 
                xlab("Number of Observations") + ylab("Change-point") + 
                scale_x_continuous(breaks = seq(1, n, floor(n/10))) + 
                scale_y_continuous(breaks = c(0, 1))
        }
        else {
            plotTsum + geom_point() + theme_bw(base_size = 16) + 
                xlab("Number of Observations") + ylab("Change-point") + 
                scale_x_continuous(breaks = breaks) + scale_y_continuous(breaks = c(0, 
                1))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
