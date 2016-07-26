\name{AIC_jgl}
\alias{AIC_jgl}
\title{
Compute Akaike Information Criterion starting from a jgl object
}
\description{
Compute Akaike Information Criterion starting from a jgl object, a vector of sample sizes and a list of covariance matrices
}
\usage{
AIC_jgl(jgl, n, S)
}
\arguments{
  \item{jgl}{
The output of \code{\link{JGL}}
}
  \item{n}{
a vector of sample sizes.
}
  \item{S}{
A list of sample covariance matrices
}
}
\details{
It uses the formula in Danaher et al. (2014)
}
\value{
An integer, the Akaike Information Criterion
}
\references{
Danaher, P., Wang, P., and Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373-397. http://doi.org/10.1111/rssb.12033
}
\author{
Giulio Costantini
}

\seealso{
\code{\link{JGL}}
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
# standardize data data within classes
sp <- split(dat[, !names(dat) == splt], dat[, splt])
sp_sc <- lapply(sp, scale)
dat <- lapply(sp_sc, data.frame)
jgl <- JGL(Y = dat, lambda1 = .1, lambda2 = .1)
S <- lapply(dat, cov)
n <- sapply(dat, nrow)
AIC_jgl(jgl, n, S)
}
}
\keyword{Joint Graphical Lasso}
\keyword{Akaike Information Criterion}
