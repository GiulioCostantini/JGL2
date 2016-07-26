# 
# # a toy example
# N <- 1000 # sample size
# sigma1 <- matrix(c(1, .5, 0, 0,
#                   .5, 1, .2, 0,
#                   0, .2, 1, 0,
#                   0, 0, 0, 1), ncol = 4)
# 
# sigma2 <- matrix(c(1, .5, .4, .4,
#                    .5, 1, .2, 0,
#                    .4, .2, 1, 0,
#                    .4, 0, 0, 1), ncol = 4)
# 
# 
# dat <- list()
# dat[[1]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma1)), Sigma = sigma1)
# dat[[2]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma2)), Sigma = sigma2)
# lapply(dat, function(x) corpcor::cor2pcor(cor(x)))
# dat <- data.frame(rbind(dat[[1]], dat[[2]]))
# dat$splt <- c(rep(1, N), rep(2, N))
# 
# aics <- JGL_AIC_surfaceplot(dat, splt = "splt", l1max = .01, l2max = .01, ncand = 10, ncores = 7, plot = TRUE)
# scatter3d(aic ~ lambda1 + lambda2, data = aics, fit = "smooth", surface = FALSE)
# 

# compute the criterioin defined by aicfun on ncand l1 x ncand l2 pairs of l1 and l2, optionally does a 3d plot of the criterion
# on the l1 and l2 values
JGL_AIC_surfaceplot <- function(dat, splt, l1min = 0, l1max = 1, l2min = 0, l2max = 1, ncand = 20, aicfun = AIC_jgl, plot = TRUE, ncores = 1, ...)
  
  # dat = a dataframe
  # splt = the column of the dataframe dat that defines multiple classes
  # l1min, l1max, l2min, l2max = minimum/maximum values for lambdas
  # ncand = number of values for lambda 1 and for lambda 2
  # aicfun = AIC / BIC / others that take as input (1) the output of JGL, (2) n = vector of sample sizes, (3) S = list of covariance matrices and give in output a criterion to minimize
  # plot = if TRUE make a 3d surface plot of the aicfun values over the l1 and l2 parameters
  # ncores = how many cores to use.
  # ... = Other parameters to be passed to JGL
{
  # standardize data data within classes
  sp <- split(dat[, !names(dat) == splt], dat[, splt])
  sp_sc <- lapply(sp, scale)
  dat <- lapply(sp_sc, data.frame)
  
  S <- lapply(dat, cov)
  n <- sapply(dat, nrow)
  
  l1cand <- seq(from = l1min, to = l1max, length.out = ncand)
  l2cand <- seq(from = l2min, to = l2max, length.out = ncand)

  lambdas <- data.frame(expand.grid(l1cand, l2cand))
  fun <- function(lambda, ...)
  {
    jgl <- JGL(Y = dat, lambda1 = lambda[1], lambda2 = lambda[2], return.whole.theta = TRUE, ...)
    aicfun(jgl = jgl, n = n, S = S)
  }
  
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun", "lambdas", "dat", "JGL", "aicfun", "n", "S"), envir = environment())
    aics <- parApply(cl = cl, lambdas, 1, fun)  
    stopCluster(cl)
  } else {
    aics <- apply(lambdas, 1, fun)  
  }

  aics <- cbind(lambdas, aics)
  names(aics) <- c("lambda1", "lambda2", "aic")

  if(plot) scatter3d(aic ~ lambda1 + lambda2, data = aics, fit = "smooth")
  
  aics
}