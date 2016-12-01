 # a toy example
# N <- 1000 # sample size
# sigma1 <- matrix(c(1, .5, 0, 0,
#                    .5, 1, .2, 0,
#                    0, .2, 1, 0,
#                    0, 0, 0, 1), ncol = 4)
# 
# sigma2 <- matrix(c(1, .5, .4, .4,
#                    .5, 1, .2, 0,
#                    .4, .2, 1, 0,
#                    .4, 0, 0, 1), ncol = 4)
# 
# dat <- list()
# dat[[1]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma1)), Sigma = sigma1)
# dat[[2]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma2)), Sigma = sigma2)
# lapply(dat, function(x) corpcor::cor2pcor(cor(x)))
# dat <- data.frame(rbind(dat[[1]], dat[[2]]))
# dat$splt <- c(rep(1, N), rep(2, N))
# system.time(
# cg1 <- JGL_AIC_widesearch(dat = dat, splt = "splt", ncand = 10, ncores = 7, optmethod = "CG")
# )

# This function performs a grid search over ncand values of l1 and l2. For each pair of values, it tries to find the minimum value of a criteiron (AIC/BIC)
# using optimization. Since the problem seems not convex, considering many starting values for l1 and l2 helps escaping local minima,
# however it is also computationally *very* intensive
JGL_AIC_widesearch <- function(dat, splt, return.whole.theta = TRUE, l1min = 0, l1max = 1, l2min = 0, l2max = 1, ncand = 20, criterion = c("ebic","aic"), gamma = 0.5, optmethod = "CG", ncores = 1, ...)
  
  # dat = a dataframe
  # splt = the column of the dataframe dat that defines multiple classes
  # return.whole.theta = see JGL
  # l1min, l1max, l2min, l2max = minimum/maximum values for lambdas
  # ncand = number of candidate lambdas initially considered (ncand l1 x ncand l2 will be considered)
  # aicfun = AIC / BIC / others that take as input (1) the output of JGL, (2) n = vector of sample sizes, (3) S = list of covariance matrices and give in output a criterion to minimize
  # optmethod = parameter method to be passed to optim
  # ncores = the function is now parallelized for windows systems. Ncores determines how many cores are used. I don't know hot it performs on other OSs, e.g., linux, mac.
  # ... = Other parameters to be passed to JGL (may not work well for parallelized code!!)
{
  
  criterion = match.arg(criterion)
  aicfun = switch(criterion, aic = AIC_jgl, ebic = BIC_jgl)
  

  # standardize data data within classes
  sp <- split(dat[, !names(dat) == splt], dat[, splt])
  sp_sc <- lapply(sp, scale)
  dat <- lapply(sp_sc, data.frame)
  
  S <- lapply(dat, cov)
  n <- sapply(dat, nrow)
  
  ###############################################
  # define the candidate pairs of lambda values #
  ###############################################
  maxl1 <- function(S) max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  l1_theormax  <- sapply(S, maxl1)
  l1max <- min(l1_theormax, l1max)
  if(l1max < l1min) stop("the specified minimum value for l1 is less than the maximun useful value for l1")
  
  # define an equally spaced sequence of ncand candidates, equally spaced between l1min and l1max
  l1cand <- seq(from = l1min, to = l1max, length.out = ncand)
  
  # Across the candidate l1, find the lowest value of lambda2 that, given that l1, makes all the elements of the precision matix
  #  equal among different classes (max difference allowed = tol), using the stats::optimize function 
  fun <- function(lambda2, lambda1, tol = 1e-5)
  {
    maxdiff <- c()
    jgl <- JGL(Y = dat, lambda2 = lambda2, lambda1 = lambda1, return.whole.theta = TRUE)
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
  
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun", "l2max", "l1cand", "JGL", "dat"), envir = environment())
    l2_theormax <- parSapply(cl = cl, 1:ncand, function(i) optimize(f = fun, lambda1 = l1cand[[i]], interval = c(0, l2max))$minimum)
    stopCluster(cl)
  } else {
    l2_theormax <-  sapply(1:ncand, function(i) optimize(f = fun, lambda1 = l1cand[[i]], interval = c(0, l2max))$minimum)
  } 
  
  # for each lambda1, define a sequence of ncand candidates, equally spaced between l2min and the updated l2max
  l2s <- sapply(l2_theormax, function(x) seq(from = l1min, to = x, length.out = ncand))
  l2s <- sapply(l2s, function(x) x)
  
  # candidate pairs of lambdas
  lambdas <-cbind(rep(l1cand, each = ncand), l2s)
  
  ######################
  # optimization stage #
  ######################
  
  # use the pairs of candidate lambda values as starting values for the optimizer
  fun <- function(lambda, ...)
  {
    averylargenumber <- 1e10
    if(any(lambda < 0)) return(averylargenumber)
    jgl <- JGL(Y = dat, lambda1 = lambda[1], lambda2 = lambda[2], return.whole.theta = TRUE, ...)
    aicfun(jgl = jgl, n = n, S = S, gamma = gamma)
  }
  
  if(ncores > 1)
  {
    cl <- makeCluster(ncores)
    clusterExport(cl, list("fun", "lambdas", "aicfun", "JGL", "dat", "n", "S", "gamma"), envir = environment())
    opt <- parApply(cl = cl, lambdas, 1, optim, fun, method = optmethod)
    stopCluster(cl)
  }  else {
    opt <- apply(lambdas, 1, optim, fun, method = optmethod)  
  }

  crits <- data.frame(t(sapply(opt, function(x) c(x$par, x$value))))
  names(crits) <- c("lambda1", "lambda2", "aic")
  
  best <- which.min(sapply(opt, function(x) x$value))
  opt <- opt[[best]]
  lambda1 <- opt$par[1]
  lambda2 <- opt$par[2]
  
  # compute the networks using the parameters lambda 1 and lambda2 choosen according to AIC
  jgl <- JGL(Y = dat, lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = return.whole.theta, ...) 
  names(jgl$theta) <- names(dat)
  
  aic <- aicfun(jgl, n, S, gamma = gamma)
  
  # return 
  list("jgl" = jgl,
       "lambda1" = lambda1,
       "lambda2" = lambda2,
       "l1theormax" = l1_theormax,
       "l2theormax" = l2_theormax,
       "aic" = aic,
       "critsFor3dPlot" = crits)
}


