\name{ComputeCost}
\alias{ComputeCost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute cost matrix for each interval
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a part of DP approach.
This function computes cost matrix for each interval.
}
\usage{
ComputeCost(x, y, n, p)
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
total number of samples
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data
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
function (x, y, n, p) 
{
    c <- matrix(0, n, n)
    for (i in 1:n) {
        A <- matrix(0, p, p)
        B <- matrix(0, p, 1)
        for (j in i:n) {
            A <- A + x[j, ] \%*\% t(x[j, ])
            B <- B + y[j] * x[j, ]
            alpha <- ginv(A) \%*\% B
            tmp <- 0
            for (t in i:j) {
                tmp <- tmp + (y[t] - (t(x[t, ]) \%*\% alpha)[1])^2
            }
            c[i, j] <- tmp
        }
    }
    return(c)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
