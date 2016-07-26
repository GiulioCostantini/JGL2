\name{JGL_cv}
\alias{JGL_cv}
\title{
Joint Graphical Lasso with k-fold crossvalidation
}
\description{
Fits a Joint Graphical Lasso model. The lambda 1 and lambda 2 parameters are identified through k-fold crossvalidation. The cost function for the k-fold crossvalidation procedure is the average predictive negative loglikelihood, as defined in Guo et al. (2011, p.5)
}
\usage{
JGL_cv(dat, splt, ncand = 20, l2max = 10, seed = 1, k = 10, ncores = 1, aicfun = AIC_jgl, ...)
}
\arguments{
\item{dat}{
A dataset that includes the variables on which the gaussian graphical models should be computed, plus an additional factor \code{"splt"} which defines different classes.
}
  \item{splt}{
Character string. The name of the variable in dat that defines different classes
}
\item{ncand}{
Integer. number of values for lambda 1 and for lambda 2
}
\item{l2max}{
Numeric. Maximum value of the lasso parameter lambda2
}
  \item{seed}{
seed parameter that ensures the exact reproducibility of the results obtained through the k-fold crossvalidation 
}
  \item{k}{
number of splits for the k-fold crossvalidation
}
  \item{ncores}{
Number of cores to use. The function is optimized for parallel computing in Windows, parallel computing may not work on other systems.
}
\item{aicfun}{
\code{AIC.jgl} or another function that take as input (1) the output of \code{JGL}, (2) n = vector of sample sizes, (3) S = list of covariance matrices. Differnetly from \code{\link{JGL_AIC_sequentialsearch}} and from \code{\link{JGL_AIC_widesearch}}, the aicfun is not used in the choice of the parameters lambda 1 and lambda 2, it has been included here only for consistency with the other functions.
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
}
\references{
Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373-397. http://doi.org/10.1111/rssb.12033
Guo, J., Levina, E., Michailidis, G., & Zhu, J. (2011). Joint estimation of multiple graphical models. Biometrika, 98(1), 1-15. http://doi.org/10.1093/biomet/asq060
}
\author{
Giulio Costantini
}
\seealso{
\code{\link{JGL}}, \code{\link{JGL_AIC_sequentialsearch}}, \code{\link{JGL_AIC_widesearch}}, \code{\link{JGL_AIC_surfaceplot}}
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
splt <- "splt"

x <- JGL_cv(dat = dat, splt = splt, ncand = 20, l2max = 1, seed = 1, k = 10, ncores = 7, aicfun = AIC)
x
}
}

\keyword{Joint Graphical Lasso}
\keyword{Crossvalidation}
