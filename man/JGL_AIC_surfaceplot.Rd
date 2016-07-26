\name{JGL_AIC_surfaceplot}
\alias{JGL_AIC_surfaceplot}
\title{
3D plot of AIC values
}
\description{
Compute the criterioin defined by \code{aicfun} on ncand l1 x ncand l2 pairs of l1 and l2. Optionally does a 3d plot of the criterion on the l1 and l2 values

}
\usage{
JGL_AIC_surfaceplot(dat, splt, l1min = 0, l1max = 1, l2min = 0, l2max = 1, ncand = 20, aicfun = AIC_jgl, plot = TRUE, ncores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{
A dataset that includes the variables on which the gaussian graphical models should be computed, plus an additional factor \code{"splt"} which defines different classes.
}
  \item{splt}{
Character string. The name of the variable in dat that defines different classes
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
  \item{aicfun}{
\code{AIC_jgl} or another function that take as input (1) the output of \code{JGL}, (2) n = vector of sample sizes, (3) S = list of covariance matrices and give in output a criterion to minimize
}
  \item{plot}{
if \code{TRUE} make a 3d surface plot of the aicfun values over the l1 and l2 parameters, using function \code{\link{scatter3d}}
}
  \item{ncores}{
Number of cores to use. The function is optimized for parallel computing in Windows, parallel computing may not work on other systems.
}
  \item{\dots}{
Additional parameters for \code{\link{JGL}}
}
}

\value{
A dataframe with the following objects
\item{lambda1}{Values of lambda 1 considered}
\item{lambda2}{Values of lambda 2 considered}
\item{aic}{The AIC or other criterion for each pair of lambda 1 and lambda 2 values}
}
\author{
Giulio Costantini
}
\examples{
\dontrun{
# # a toy example
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

aics <- JGL_AIC_surfaceplot(dat, splt = "splt", l1max = .01, l2max = .01, ncand = 10, ncores = 7, plot = TRUE)
scatter3d(aic ~ lambda1 + lambda2, data = aics, fit = "smooth", surface = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Joint Graphical Lasso}% use one of  RShowDoc("KEYWORDS")
\keyword{3D}% __ONLY ONE__ keyword per line
