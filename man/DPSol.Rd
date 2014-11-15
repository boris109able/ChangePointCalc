\name{DPSol}
\alias{DPSol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute change-points using DP approach.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a part of DP approach.
This fucntion computes change-points using DP approach.
}
\usage{
DPSol(x, y, n, p, Kstar, c)
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
total number of data
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data
}
  \item{Kstar}{
%%     ~~Describe \code{Kstar} here~~
given number of change-points
}
  \item{c}{
%%     ~~Describe \code{c} here~~
cost matrix c
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
Bingwen
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
function (x, y, n, p, Kstar, c) 
{
    result <- DPCalc(x, y, n, p, Kstar, c)
    r <- result$r
    s <- result$s
    K <- Kstar + 1
    sol <- vector()
    sol <- append(sol, n)
    while (K > 0) {
        sol <- append(sol, s[K, n])
        n <- s[K, n]
        K <- K - 1
    }
    sol2 <- vector()
    len <- length(sol)
    for (i in len:1) {
        sol2 <- append(sol2, sol[i])
    }
    return(sol2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
