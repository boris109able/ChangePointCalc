\name{DPmain}
\alias{DPmain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Main function of DP approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Main function of DP approach
}
\usage{
DPmain(x, y, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
input data x
}
  \item{y}{
%%     ~~Describe \code{y} here~~
input data y
}
  \item{K}{
%%     ~~Describe \code{K} here~~
given number of change-points
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
function (x, y, K) 
{
    n <- dim(x)[1]
    p <- dim(x)[2]
    c1 <- computeC(x, y, n, p)
    sol <- DPSol(x, y, n, p, K + 1, c1)
    predictedAlpha <- PrintAlpha(x, y, sol)
    return(list(Sol = sol, AlphaHat <- predictedAlpha))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
