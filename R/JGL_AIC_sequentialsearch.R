# a toy example
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
# JGL_AIC_sequentialsearch(dat = dat, splt = "splt", ncores = 1)

# parallelization does not seem doing good here, perhaps the example is too simple!
# system.time(
#  JGL_AIC_sequentialsearch(dat = dat, splt = "splt", ncores = 7)
# )
# 
# system.time(
#  JGL_AIC_sequentialsearch(dat = dat, splt = "splt", ncores = 1)
# )

# performs a sequential search of the values of l1 and l2. First it finds the best fitting value of l1 according to the selected aicfun (AIC/BIC)
# then it finds the best value of l2 given the l1 selected at the first stage.
# Optionally, it can use the values of l1 and l2 as the input for a 2-dimensional optimization, to search if the aicfun can be further improved.

JGL_AIC_sequentialsearch <- function(dat, splt, return.whole.theta = TRUE, l1min = 0, l1max = 1, l2min = 0, l2max = 1, ncand = 20, aicfun = AIC_jgl, globalopt = TRUE, optmethod = "CG", ncores = 1, ...)
  
  # dat = a dataframe
  # splt = the column of the dataframe dat that defines multiple classes
  # return.whole.theta = see JGL
  # l1min, l1max, l2min, l2max = minimum/maximum values for lambdas
  # ncand = number of candidate lambdas initially considered
  # aicfun = AIC / BIC / others that take as input (1) the output of JGL, (2) n = vector of sample sizes, (3) S = list of covariance matrices and give in output a criterion to minimize
  # globalopt = if TRUE, perform a further step that tries a 2-dimensional optimization of both l1 and l2 using optim optimizer
  # optmethod = parameter method to be passed to optim
  # ... = Other parameters to be passed to JGL (may not work well for parallelized code!!)
{
  # standardize data data within classes
  sp <- split(dat[, !names(dat) == splt], dat[, splt])
  sp_sc <- lapply(sp, scale)
  dat <- lapply(sp_sc, data.frame)

  S <- lapply(dat, cov)
  n <- sapply(dat, nrow)
  
  ##################
  # Select lambda1 #
  ##################
  
  # Define a range of candidate values for lambda 1. Select, as the largest lambda1 (l1max)
  # the one that would make all edges = 0 in at least one of the networks if the glasso
  # penalty was applied alone (i.e., with lambda2 = 0). 
  # code for l1max is adapted from package qgraph, which in turn took it from package huge
  maxl1 <- function(S) max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  l1_theormax  <- sapply(S, maxl1)
  l1max <- min(l1_theormax, l1max)
  if(l1max < l1min) stop("the specified minimum value for l1 is less than the maximun useful value for l1")
  
  # define an  equally spaced sequence of ncand candidates
  l1cand <- seq(from = l1min, to = l1max, length.out = ncand)
  
  # select the value of l1 with minimal AIC among the candidates, considering a minimal value of l2
  fun <- function(lambda1, ...)
  {
    jgl <- JGL(Y = dat, lambda1 = lambda1, lambda2 = l2min, return.whole.theta = TRUE, ...)
    aicfun(jgl = jgl, n = n, S = S)
  }
  
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun", "dat", "n", "S", "JGL", "l1cand", "l2min", "aicfun"), envir = environment())
    aicsl1 <- parSapply(cl = cl, l1cand, fun)  
    stopCluster(cl)
  } else {
    aicsl1 <- sapply(l1cand, fun)  
  }
  
  l1cand <- c(l1min, l1cand, l1max)
  l1min <- l1cand[which.min(aicsl1)]
  l1max <- l1cand[which.min(aicsl1) + 2]
  # see if AIC can be further improved in the neighborhood of the chosen value
  opt <- optimize(f = fun, interval = c(l1min, l1max))
  lambda1 <- opt$minimum
  
  ##################
  # Select lambda2 #
  ##################
  
  #  find the lowest value of lambda2 that, given l1, that makes all the elements of the precision matix
  #  equal among different classes (max difference allowed = tol), using the stats::optimize function 
  
  fun <- function(lambda2, tol = 1e-5, ...)
  {
    maxdiff <- c()
    jgl <- JGL(Y = dat, lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = TRUE, ...)
    thetas <- jgl$theta
    for(i in 2:length(thetas))
    {
      for(j in 1:(i-1))
      {
        maxdiff <- c(maxdiff, max(abs(thetas[[i]] - thetas[[j]])))
      }
    }
    maxdiff <- max(maxdiff)
    
    ifelse(maxdiff > tol, l2max+1, lambda2)
  }
  
  opt <- optimize(f = fun, interval = c(0, l2max)) 
  
  l2_theormax <- opt$minimum
  l2max <- min(c(l2max, l2_theormax))
  
  # define a sequence of ncand candidate values for l2, equally spaced between l2min and l2max
  l2cand <- seq(from = l2min, to = l2max, length.out = ncand)
  
  # select the value of l2 with minimal AIC among the candidates
  fun <- function(lambda2, ...)
  {
    jgl <- JGL(Y = dat, lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = TRUE, ...)
    aicfun(jgl = jgl, n = n, S = S)
  }
  
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun", "dat", "lambda1", "n", "S", "JGL", "l2cand", "aicfun"), envir = environment())
    aicsl2 <- parSapply(cl = cl, l2cand, fun) 
    stopCluster(cl)
  } else {
    aicsl2 <- sapply(l2cand, fun)
  }
  
  # find the lambda 2 that locally minimizes AIC around the choosen value
  l2cand <- c(l2min, l2cand, l2max)
  l2min <- l2cand[which.min(aicsl2)]
  l2max <- l2cand[which.min(aicsl2) + 2]
  # see if AIC can be further improved in the neighborhood of the chosen value
  opt <- optimize(f = fun, interval = c(l2min, l2max))
  lambda2 <- opt$minimum
  
  if(globalopt)
  {
    startl1 <- lambda1
    startl2 <- lambda2
    averylargenumber <- 1e10
    fun <- function(lambda, ...)
    {
      if(any(lambda < 0)) return(averylargenumber)
      jgl <- JGL(Y = dat, lambda1 = lambda[1], lambda2 = lambda[2], return.whole.theta = TRUE, ...)
      aicfun(jgl = jgl, n = n, S = S)
    }
    
    opt <- optim(par = c(startl1, startl2), fun, method = "BFGS")
    
    lambda1 <- opt$par[1]
    lambda2 <- opt$par[2]
  }
  
   # compute the networks using the parameters lambda 1 and lambda2 choosen according to AIC
  jgl <- JGL(Y = dat, lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = return.whole.theta, ...)
  
  names(jgl$theta) <- names(dat)
  aic <- aicfun(jgl, n, S)
  
  # return 
  list("jgl" = jgl,
       "lambda1" = lambda1,
       "lambda2" = lambda2,
       "l1theormax" = l1_theormax,
       "l2theormax" = l2_theormax,
       "aic" = aic)
}

