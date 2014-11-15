\name{oneDim}
\alias{oneDim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The function to calculate SGL of linear regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To implement interface to calculate SGL using C++
}
\usage{
oneDim(data, index, thresh = 0.01, nlam = 20, lambdas = NULL, beta.naught = rep(0, ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 0.01, gamma = 0.8, step = 1, reset = 10, alpha = 0.95, min.frac = 0.05, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
list of input data x and y
}
  \item{index}{
%%     ~~Describe \code{index} here~~
group index
}
  \item{thresh}{
%%     ~~Describe \code{thresh} here~~
threshold 
}
  \item{nlam}{
%%     ~~Describe \code{nlam} here~~
number of lambdas long the regularization path
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
list of values of lambdas set by user
}
  \item{beta.naught}{
%%     ~~Describe \code{beta.naught} here~~
beta.naught
}
  \item{inner.iter}{
%%     ~~Describe \code{inner.iter} here~~
innter.iter
}
  \item{outer.iter}{
%%     ~~Describe \code{outer.iter} here~~
outer.iter
}
  \item{outer.thresh}{
%%     ~~Describe \code{outer.thresh} here~~
threshold of outer loop
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
gamma
}
  \item{step}{
%%     ~~Describe \code{step} here~~
step
}
  \item{reset}{
%%     ~~Describe \code{reset} here~~
reset
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
alpha
}
  \item{min.frac}{
%%     ~~Describe \code{min.frac} here~~
min.frac
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
verbose
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
Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011) A Sparse-Group Lasso,
http://www-stat.stanford.edu/~nsimon/SGL.pdf
}
\author{
%%  ~~who you are~~
Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani
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
function (data, index, thresh = 0.01, nlam = 20, lambdas = NULL, 
    beta.naught = rep(0, ncol(data$x)), inner.iter = 100, outer.iter = 100, 
    outer.thresh = 0.01, gamma = 0.8, step = 1, reset = 10, alpha = 0.95, 
    min.frac = 0.05, verbose = FALSE) 
{
    if (is.null(lambdas)) {
        lambdas <- betterPathCalc(data = data, index = index, 
            alpha = alpha, min.frac = min.frac, nlam = nlam, 
            type = "linear")
    }
    dyn.load("F://R files//PackageTest03//PackageTest//src//PackageTest.dll")
    X <- data$x
    y <- data$y
    n <- nrow(X)
    p <- ncol(X)
    ord <- order(index)
    index <- index[ord]
    X <- X[, ord]
    unOrd <- match(1:length(ord), ord)
    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0, (num.groups + 1))
    for (i in 1:num.groups) {
        range.group.ind[i] <- min(which(index == groups[i])) - 
            1
    }
    range.group.ind[num.groups + 1] <- ncol(X)
    group.length <- diff(range.group.ind)
    nlam = length(lambdas)
    beta.old <- rep(0, ncol(X))
    beta.is.zero <- rep(1, num.groups)
    beta <- array(0, c(ncol(X), nlam))
    eta <- rep(0, n)
    for (k in 1:nlam) {
        beta.is.zero <- rep(1, num.groups)
        beta.old <- rep(0, ncol(X))
        eta <- rep(0, n)
        junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), 
            index = as.integer(index), nrow = as.integer(nrow(X)), 
            ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), 
            rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), 
            lambda1 = as.double(lambdas[k] * alpha), lambda2 = as.double(lambdas[k] * 
                (1 - alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), 
            outerIter = as.integer(outer.iter), thresh = as.double(thresh), 
            outerThresh = as.double(outer.thresh), eta = as.double(eta), 
            gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), 
            step = as.double(step), reset = as.integer(reset), 
            PACKAGE = "PackageTest")
        beta.new <- junk$beta
        beta[, k] <- beta.new
        beta.is.zero <- junk$betaIsZero
        eta <- junk$eta
        beta.old <- beta.new
        if (verbose == TRUE) {
            write(paste("***Lambda", k, "***"), "")
        }
    }
    return(list(beta = beta[unOrd, ], lambdas = lambdas))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
