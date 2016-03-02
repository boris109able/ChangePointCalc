\name{SGLPlotNorm}
\alias{SGLPlotNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot change-points using norm by SGL approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot change-points using norm by SGL approach
}
\usage{
SGLPlotL1Norm(fit, n, p, num = 1, breaks = NULL, norm = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
result of SGL approach
}
  \item{n}{
%%     ~~Describe \code{n} here~~
total number of data
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data
}
  \item{num}{
%%     ~~Describe \code{num} here~~
the num'th lambda used by user
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
list of breaks where to put a tick on the x-axis
}
 \item{norm}{
%%     ~~Describe \code{breaks} here~~
norm=2 is l2 norm, norm=1 is l1 norm
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
function(fit, n, p, num=1, breaks=NULL,norm=2)
{
  # this function plot the l1 norm of theta in the result
  # num is the num th lambda we used
  # n is the number of data sample
  # p is the number of dimension
  #norm=2 is l2 norm, norm=1 is l1 norm
  #breaks is the break tick the user want to show on the plot, 
  #if breaks is NULL the default setting is to divide x-axis into 10 equal segment and label the 
  #tick accordingly
  t_sum = matrix(NA, n, 1)
  if(norm==1)
  {
    for (i in 1:n)
    { t_sum[i] = sum(abs(fit$beta[((i-1)*p+1):(i*p), num]))} 
  }
  if(norm==2)
  {
    for (i in 1:n)
    { t_sum[i] = sqrt(sum((fit$beta[((i-1)*p+1):(i*p), num])^2))}
  }
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
