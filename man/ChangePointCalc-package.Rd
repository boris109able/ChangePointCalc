\name{ChangePointCalc-package}
\alias{ChangePointCalc-package}
\alias{ChangePointCalc}
\docType{package}
\title{
Offline change-Points estimation via DP and SGL
}
\description{
Calculate change-points by DP approach and SGL approach in linear regression models.
}
\details{
\tabular{ll}{
Package: \tab ChangePointCalc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-25\cr
License: \tab GPL\cr
}

}
\author{
Bingwen Zhang
Maintainer: Bingwen Zhang <tszhang109@gmail.com>
}
\references{
Bingwen Zhang, Jun Geng and Lifeng Lai, Detecting changes in regression models via sparse group lasso,
http://users.wpi.edu/~bzhang/icassp_v1.pdf

Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011) A Sparse-Group Lasso,
http://www-stat.stanford.edu/~nsimon/SGL.pdf
}

\keyword{ package }

\examples{
ChangePoints <- c(101,201)
n <- 300
p <- 20
SData <- SynthesizeData(ChangePoints, sigma1=1, sigma2=0.1, n, p, coefficients=NULL, type="SGL")
np <- n*p
index <- ceiling(1:np/p)
data <- list(x=SData$x,y=SData$y)
resultSGL <- SGL(data,index,type='linear',alpha=0.20)
SGLPlotChangePoint(resultSGL, n, p, num=20, breaks=NULL)

SData <- SynthesizeData(ChangePoints, sigma1=1, sigma2=0.1, n, p, coefficients=NULL, type="DP")
resultDP <- DPmain(x=SData$x, y=SData$y, 3)
DPPlotChangePoint(resultDP$Sol, n, p)
}
