\name{SynthesizeData}
\alias{SynthesizeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate synthesized data used by SGL approach
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate synthesized data used by SGL approach
}
\usage{
SynthesizeData(ChangePoints, sigma1, sigma2, n = 1000, p = 20, coefficients = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ChangePoints}{
%%     ~~Describe \code{ChangePoints} here~~
list of change-points set by user
}
  \item{sigma1}{
%%     ~~Describe \code{sigma1} here~~
variance of signal
}
  \item{sigma2}{
%%     ~~Describe \code{sigma2} here~~
variance of noise
}
  \item{n}{
%%     ~~Describe \code{n} here~~
total number of data generated
}
  \item{p}{
%%     ~~Describe \code{p} here~~
dimension of data generated
}
  \item{coefficients}{
%%     ~~Describe \code{coefficients} here~~
matrix of coefficients of each interval set by user
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
function (ChangePoints, sigma1, sigma2, n = 1000, p = 20, coefficients = NULL) 
{
    np <- n * p
    KReal <- length(ChangePoints)
    TmpPoints <- c(1, ChangePoints, n + 1)
    x <- matrix(rnorm(n * p) * sigma1, ncol = p, nrow = n)
    XTilde <- matrix(NA, n, np)
    for (i in 1:n) {
        for (j in 1:i) {
            XTilde[i, ((j - 1) * p + 1):(j * p)] <- x[i, ]
        }
    }
    if (is.null(coefficients)) {
        coefficients <- matrix(NA, KReal + 1, p)
        for (i in 1:(KReal + 1)) {
            if (i\%\%2 == 1) {
                coefficients[i, ] <- 2 * rep.int(1, p)
            }
            else {
                coefficients[i, ] <- -2 * rep.int(1, p)
            }
        }
    }
    beta <- matrix(NA, n, p)
    for (i in 1:(KReal + 1)) {
        beta[TmpPoints[i]:(TmpPoints[i + 1] - 1), ] <- t(matrix(rep(coefficients[i, 
            ], TmpPoints[i + 1] - TmpPoints[i]), p, TmpPoints[i + 
            1] - TmpPoints[i]))
    }
    y <- matrix(NA, n, 1)
    for (i in 1:n) {
        y[i] <- (x[i, ] \%*\% beta[i, ])[1] + matrix(rnorm(1) * 
            sigma2, ncol = 1, nrow = 1)
    }
    data <- list(x = XTilde, y = y, coefficients = coefficients, 
        beta = beta)
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
