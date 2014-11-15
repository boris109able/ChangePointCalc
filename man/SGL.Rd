\name{SGL}
\alias{SGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
main function of SGL
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
main function of SGL
}
\usage{
SGL(data, index, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1, nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
input data, x and y
}
  \item{index}{
%%     ~~Describe \code{index} here~~
index
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
maxit
}
  \item{thresh}{
%%     ~~Describe \code{thresh} here~~
threshold
}
  \item{min.frac}{
%%     ~~Describe \code{min.frac} here~~
minfrac
}
  \item{nlam}{
%%     ~~Describe \code{nlam} here~~
number of lambdas
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
gamma
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
flag of standardization
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
verbose
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
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
sequence of lambdas
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
result of SGL
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
function (data, index, maxit = 1000, thresh = 0.001, min.frac = 0.1, 
    nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, 
    step = 1, reset = 10, alpha = 0.95, lambdas = NULL) 
{
    X.transform <- NULL
    if (standardize == TRUE) {
        X <- data$x
        means <- apply(X, 2, mean)
        X <- t(t(X) - means)
        var <- apply(X, 2, function(x) (sqrt(sum(x^2))))
        X <- t(t(X)/var)
        data$x <- X
        X.transform <- list(X.scale = var, X.means = means)
    }
    if (standardize == TRUE) {
        intercept <- mean(data$y)
        data$y <- data$y - intercept
    }
    Sol <- oneDim(data, index, thresh, inner.iter = maxit, outer.iter = maxit, 
        outer.thresh = thresh, min.frac = min.frac, nlam = nlam, 
        lambdas = lambdas, gamma = gamma, verbose = verbose, 
        step = step, reset = reset, alpha = alpha)
    if (standardize == TRUE) {
        Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, 
            intercept = intercept, X.transform = X.transform)
    }
    if (standardize == FALSE) {
        Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, 
            X.transform = X.transform)
    }
    class(Sol) = "SGL"
    return(Sol)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
