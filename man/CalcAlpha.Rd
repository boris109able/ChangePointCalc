\name{CalcAlpha}
\alias{CalcAlpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate coefficients of each interval
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a part of DP approach.
This function calculates coefficients of each interval.
}
\usage{
CalcAlpha(x, y, n, p, i, j)
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
  \item{n}{
%%     ~~Describe \code{n} here~~
number of data samples
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data
}
  \item{i}{
%%     ~~Describe \code{i} here~~
from index i
}
  \item{j}{
%%     ~~Describe \code{j} here~~
to index j
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
function (x, y, n, p, i, j) 
{
    A <- matrix(0, p, p)
    B <- matrix(0, p, 1)
    for (t in i:j) {
        A <- A + x[t, ] \%*\% t(x[t, ])
        B <- B + y[t] * x[t, ]
    }
    alpha <- ginv(A) \%*\% B
    return(alpha)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
