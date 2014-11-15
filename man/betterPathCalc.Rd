\name{betterPathCalc}
\alias{betterPathCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To calculate the path of lambdas
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is to calculate a sequence of lambdas which consist of the path.
The sequence of lambdas is returned.
}
\usage{
betterPathCalc(data, index, alpha = 0.95, min.frac = 0.05, nlam = 20, type = "linear")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
list of input data x and y
}
  \item{index}{
%%     ~~Describe \code{index} here~~
index we use
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
coefficient to balance two penalties
}
  \item{min.frac}{
%%     ~~Describe \code{min.frac} here~~
minfrac times lambda max
}
  \item{nlam}{
%%     ~~Describe \code{nlam} here~~
number of lambdas(default 20)
}
  \item{type}{
%%     ~~Describe \code{type} here~~
"linear"
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
function (data, index, alpha = 0.95, min.frac = 0.05, nlam = 20, 
    type = "linear") 
{
    reset <- 10
    step <- 1
    gamma <- 0.8
    inner.iter <- 1000
    outer.iter <- 1000
    thresh = 10^(-2)
    outer.thresh = thresh
    n <- nrow(data$x)
    if (type == "linear") {
        X <- data$x
        resp <- data$y
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
    }
    lambda.max <- rep(0, num.groups)
    if ((alpha != 0) * (alpha != 1)) {
        for (i in 1:num.groups) {
            ind <- groups[i]
            X.fit <- X[, which(index == ind)]
            cors <- t(X.fit) \%*\% resp
            ord.cors <- sort(abs(cors), decreasing = TRUE)
            if (length(ord.cors) > 1) {
                norms <- rep(0, length(cors) - 1)
                lam <- ord.cors/alpha
                for (j in 1:(length(ord.cors) - 1)) {
                  norms[j] <- sqrt(sum((ord.cors[1:j] - ord.cors[j + 
                    1])^2))
                }
                if (norms[1] >= lam[2] * (1 - alpha) * sqrt(group.length[i])) {
                  our.cors <- ord.cors[1]
                  our.range <- c(ord.cors[2], ord.cors[1])/alpha
                }
                else {
                  if (norms[length(ord.cors) - 1] <= lam[length(ord.cors)] * 
                    (1 - alpha) * sqrt(group.length[i])) {
                    our.cors <- ord.cors
                    our.range <- c(0, ord.cors[length(ord.cors)])/alpha
                  }
                  else {
                    my.ind <- max(which(norms[-length(norms)] <= 
                      lam[2:(length(norms))] * (1 - alpha) * 
                        sqrt(group.length[i]))) + 1
                    our.cors <- ord.cors[1:my.ind]
                    our.range <- c(ord.cors[my.ind + 1], ord.cors[my.ind])/alpha
                  }
                }
                nn <- length(our.cors)
                if (alpha == 0.5) {
                  alpha = 0.500001
                }
                A.term <- nn * alpha^2 - (1 - alpha)^2 * group.length[i]
                B.term <- -2 * alpha * sum(our.cors)
                C.term <- sum(our.cors^2)
                lams <- c((-B.term + sqrt(B.term^2 - 4 * A.term * 
                  C.term))/(2 * A.term), (-B.term - sqrt(B.term^2 - 
                  4 * A.term * C.term))/(2 * A.term))
                lambda.max[i] <- min(subset(lams, lams >= our.range[1] & 
                  lams <= our.range[2]))
            }
            if (length(ord.cors) == 1) {
                lambda.max[i] <- ord.cors
            }
        }
    }
    if (alpha == 1) {
        lambda.max <- abs(t(X) \%*\% resp)
    }
    if (alpha == 0) {
        for (i in 1:num.groups) {
            ind <- groups[i]
            X.fit <- X[, which(index == ind)]
            cors <- t(X.fit) \%*\% resp
            lambda.max[i] <- sqrt(sum(cors^2))/sqrt(group.length[i])
        }
    }
    max.lam <- max(lambda.max)
    min.lam <- min.frac * max.lam
    lambdas <- exp(seq(log(max.lam), log(min.lam), (log(min.lam) - 
        log(max.lam))/(nlam - 1)))
    return(lambdas/nrow(X))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
