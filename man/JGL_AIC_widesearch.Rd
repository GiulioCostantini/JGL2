\name{JGL_AIC_widesearch}
\alias{JGL_AIC_widesearch}
\title{
Fits a Joing Graphical Lasso model. The lambda 1 and lambda 2 parameters are identified by optimizing AIC (or another criterion passed to aicfun).
}
\description{
Find the minimum value of a criteiron defined by aicfun (e.g., AIC) using 2-dimensional optimization. Since the problem is not convex, considering many starting values for l1 and l2 helps escaping local minima,
therefore ncand values of lambda1 x ncand values of lambda2 are considered. This function should perform better than \code{\link{JGL_AIC_sequentialsearch}}, however it is also much more computationally intensive.
The function is optimized for parallel computing on Windows machines. on Other OSs parallelization may not work, in this case, use ncores = 1.
}
\usage{
JGL_AIC_widesearch(dat, splt, return.whole.theta = TRUE, l1min = 0, l1max = 1, l2min = 0, l2max = 1, ncand = 20, criterion = c("ebic","aic"), gamma = 0.5, optmethod = "CG", ncores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{
A dataset that includes the variables on which the gaussian graphical models should be computed, plus an additional factor \code{"splt"} which defines different classes.
}
  \item{splt}{
Character string. The name of the variable in dat that defines different classes
}
  \item{return.whole.theta}{
Logical. Parameter passed directly to JGL, determines whether the whole concentration matrix should be returned in output. See \code{\link{JGL}}
}
  \item{l1min}{
Numeric. Minimum value of the lasso parameter lambda1
}
  \item{l1max}{
Numeric. Maximum value of the lasso parameter lambda1
}
  \item{l2min}{
Numeric. Minimum value of the lasso parameter lambda2
}
  \item{l2max}{
Numeric. Maximum value of the lasso parameter lambda2
}
  \item{ncand}{
Integer. number of values for lambda 1 and for lambda 2
}
  \item{criterion}{
The criterion for selecting the lambda values. Can be "aic" for the Akaike information criterion (AIC) or "ebic" for the Extended Bayes Information Criterion (EBIC).
}
  \item{gamma}{
The gamma value for the EBIC criterion. A value of 0 results in the BIC.
}
  \item{optmethod}{
character string. Determines the optimization method. See parameter \code{method} in function \code{\link{optim}}
}
  \item{ncores}{
Number of cores to use. The function is optimized for parallel computing in Windows, parallel computing may not work on other systems.
}
\item{\dots}{
Other parameters for \code{link{JGL}}. This feature is still in development, it may not work as intended, especially if using parallel computing.
}

}

\value{
\item{jgl}{The output of \code{\link{JGL}}}
\item{lambda1}{the value of lambda 1 that minimizes the output of aicfun}
\item{l1theormax}{the minimal value of lambda 1 that would result in at least one completely disconnected network (all missing edges). Values of lambda 1 > l1theormax are not considered}
\item{l2theormax}{for each candidate value of lambda 1, the minimal value of lambda 2 that would make the networks in the different classes all equal to each other. For each lambda 1, values of lambda 2 > l2theormax are not considered}
\item{aic}{The Akaike Information Criterion (or another criterion to minimize passed to \code{aicfun}) corresponding to the solution in \code{jgl}}
\item{critsFor3dPlot}{Please ignore this output, only for internal verification purposes. A list of the lambda 1 and 2 values identified from each starting pair, as well as the corresponding AIC}
}
\references{
Danaher, P., Wang, P., and Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373-397. http://doi.org/10.1111/rssb.12033
}
\author{
Giulio Costantini
}

\seealso{
\code{\link{JGL}}, \code{\link{JGL_AIC_sequentialsearch}}, \code{\link{JGL_cv}}, \code{\link{JGL_AIC_surfaceplot}}
}
\examples{
\dontrun{
N <- 1000 # sample size
sigma1 <- matrix(c(1, .5, 0, 0,
                   .5, 1, .2, 0,
                   0, .2, 1, 0,
                   0, 0, 0, 1), ncol = 4)

sigma2 <- matrix(c(1, .5, .4, .4,
                   .5, 1, .2, 0,
                   .4, .2, 1, 0,
                   .4, 0, 0, 1), ncol = 4)

dat <- list()
dat[[1]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma1)), Sigma = sigma1)
dat[[2]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma2)), Sigma = sigma2)
lapply(dat, function(x) corpcor::cor2pcor(cor(x)))
dat <- data.frame(rbind(dat[[1]], dat[[2]]))
dat$splt <- c(rep(1, N), rep(2, N))
system.time(
cg1 <- JGL_AIC_widesearch(dat = dat, splt = "splt", ncand = 10, ncores = 7, optmethod = "CG")
)
}
}
\keyword{Joint Graphical Lasso}
\keyword{AIC}