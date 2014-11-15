\name{DPPlotL1Norm}
\alias{DPPlotL1Norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot l1-norm of coefficients of each interval computed by DP approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot l1-norm of coefficients of each interval computed by DP approach
}
\usage{
DPPlotL1Norm(predictedAlpha, sol, n, p, breaks = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictedAlpha}{
%%     ~~Describe \code{predictedAlpha} here~~
predicted coefficients solved by DP approach
}
  \item{sol}{
%%     ~~Describe \code{sol} here~~
list of change-points computed by DP approach
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
function (predictedAlpha, sol, n, p, breaks = NULL) 
{
    predictedTheta <- matrix(0, length(sol) - 1, p)
    for (i in 1:(length(sol) - 1)) {
        if (i == 1) {
            predictedTheta[i, ] <- predictedAlpha[i, ]
        }
        else {
            predictedTheta[i, ] <- predictedAlpha[i, ] - predictedAlpha[i - 
                1, ]
        }
    }
    t_sum = mat.or.vec(n, 1)
    for (i in 1:(length(sol) - 1)) {
        t_sum[sol[i]] <- sum(abs(predictedTheta[i, ]))
    }
    plotData <- data.frame(pos = 1:length(t_sum), normL1 = t_sum)
    plotTsum <- ggplot(plotData, aes(pos, normL1))
    {
        if (is.null(breaks)) {
            plotTsum + geom_point() + theme_bw(base_size = 16) + 
                xlab("Number of Observations") + ylab(expression(paste(group("|", 
                group("|", theta[i], "|"), "|"))[1])) + scale_x_continuous(breaks = seq(1, 
                n, floor(n/10)))
        }
        else {
            plotTsum + geom_point() + theme_bw(base_size = 16) + 
                xlab("Number of Observations") + ylab(expression(paste(group("|", 
                group("|", theta[i], "|"), "|"))[1])) + scale_x_continuous(breaks = breaks)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
