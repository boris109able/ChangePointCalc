\name{PrintAlpha}
\alias{PrintAlpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
calculate coefficients of each interval
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a part of DP approach.
This function calculate coefficients of each interval.
}
\usage{
PrintAlpha(x, y, sol)
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
  \item{sol}{
%%     ~~Describe \code{sol} here~~
list of change-points by DP approach
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
function (x, y, sol) 
{
    n <- dim(x)[1]
    p <- dim(x)[2]
    K <- length(sol)
    Alpha <- matrix(NA, K - 1, p)
    if (K == 2) {
        Alpha[i, ] <- CalcAlpha(x, y, n, p, sol[1], sol[2])
    }
    else {
        for (i in 1:(K - 1)) {
            if (i == (K - 1)) {
                Alpha[i, ] <- CalcAlpha(x, y, n, p, sol[i], sol[i + 
                  1])
            }
            else {
                Alpha[i, ] <- CalcAlpha(x, y, n, p, sol[i], sol[i + 
                  1] - 1)
            }
        }
    }
    return(Alpha)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
