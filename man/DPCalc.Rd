\name{DPCalc}
\alias{DPCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate change-points using DP approach given number of change-points.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a part of DP approach.
This function calculates change-points using DP approach given number of change-points.
}
\usage{
DPCalc(x, y, n, p, Kstar, c)
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
the given number of change-points
}
  \item{c}{
%%     ~~Describe \code{c} here~~
cost matrix c (can be computed by ComputeCost() function)
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
function (x, y, n, p, Kstar, c) 
{
    r <- matrix(0, ncol = n, nrow = Kstar + 1)
    s <- matrix(0, ncol = n, nrow = Kstar + 1)
    for (t in 1:n) {
        r[1, t] <- c[1, t]
        s[1, t] <- 1
    }
    for (K in 2:(Kstar + 1)) {
        for (t in K:n) {
            q = -1
            for (j in (K - 1):(t - 1)) {
                tmp <- r[K - 1, j] + c[j + 1, t]
                if (q > tmp || q < 0) {
                  q <- tmp
                  s[K, t] <- j + 1
                }
            }
            r[K, t] <- q
        }
    }
    return(list(r = r, s = s))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
