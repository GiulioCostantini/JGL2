# # a toy example
# N <- 100 # sample size
# sigma1 <- matrix(c(1, .5, 0, 0,
#                   .5, 1, .2, 0,
#                   0, .2, 1, 0,
#                   0, 0, 0, 1), ncol = 4)
# 
# sigma2 <- matrix(c(1, .5, 0, .4,
#                    .5, 1, .2, 0,
#                    0, .2, 1, 0,
#                    .4, 0, 0, 1), ncol = 4)
# 
# 
# dat <- list()
# dat[[1]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma1)), Sigma = sigma1)
# dat[[2]] <- MASS::mvrnorm(n = N, mu = rep(0, ncol(sigma2)), Sigma = sigma2)
# lapply(dat, function(x) corpcor::cor2pcor(cor(x)))
# dat <- data.frame(rbind(dat[[1]], dat[[2]]))
# dat$splt <- c(rep(1, N), rep(2, N))
# splt <- "splt"
# 
# x <- JGL_cv(dat = dat, splt = splt, ncand = 20, l2max = 1, seed = 1, k = 10, ncores = 7, aicfun = AIC)
# x


  
# select lasso parameters through k-fold crossvalidation
JGL_cv <- function(dat, splt, ncand = 20, l2max = 10, seed = 1, k = 10, ncores = 1, aicfun = AIC_jgl, ...)
  
  # dat = a dataframe
  # splt = the column of the dataframe dat that defines multiple classes
  # l2max = maximum value considered for lambda 2, select a reasonably value
  # ncand = number of candidate lambdas initially considered (ncand l1 x ncand l2 will be considered). The computational time increases with ncand^2
  # seed  = random seed for k-fold cross validation
  # k = number of splits for the k-fold crossvalidation
  # ncores = number of pc cores to use. This function is optimized to use many cores on a windows PC, I don't know whether parallelization works on other OSs such as linux/mac.
  # ... = Other parameters to be passed to JGL (may not work well for parallelized code!!)
{
  # standardize data data within classes
  sp <- split(dat[, !names(dat) == splt], dat[, splt])
  sp_sc <- lapply(sp, scale)
  dat <- lapply(sp_sc, data.frame)
  
  S <- lapply(dat, cov)
  n <- sapply(dat, nrow)

  # find l1_theormax, the value of l1 that would cause all edges in at least one network to be missing
  maxl1 <- function(S) max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  l1_theormax  <- min(sapply(S, maxl1)) 
  
  # define an equally  spaced sequence of candidate l1, between 0 and l1_theormax. 
  l1cand <- seq(from = 0, to = l1_theormax, length.out = ncand)

  ######################################
  # find candidate values of l1 and l2 #
  ######################################

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
  
  # for each lambda1, define a sequence of ncand candidates, equally spaced between 0 and the updated l2max
  l2s <- sapply(l2_theormax, function(x) seq(from = 0, to = x, length.out = ncand))
  l2s <- sapply(l2s, function(x) x)
  
  # candidate pairs of lambdas
  lambdas <-data.frame(cbind(rep(l1cand, each = ncand), l2s))
  names(lambdas) <- c("l1", "l2")
  
  ##################################################################################
  # Perform a search of the values of lambda 1 and 2 using k-fold cross validation #
  ##################################################################################
  # some of the code for the k-fold cross validation has been taken from the R package parcor, function adalasso
  n <- sapply(dat, nrow)
  nclasses <- length(dat)
  if(!missing(seed)) set.seed(seed)
  # all 
  all.folds <- lapply(n, function(x) split(sample(1:x), rep(1:k, length = x)))
  
  crossval <- function(i, all.folds = all.folds, nclasses = nclasses, ncores = ncores, lambdas = lambdas, ...)
  {
    # for a signel fold, function crossval computes Sigmas, the covaraince matrices of the training set for each class,
    # and Omegas, a list that includes, for each pair of lambdas, the estimated concentration matrices for each class
    # finally, using function anpll, for each pair of lambdas it computes the average negative predictive log likelihood function 
    # in each training sample
    # i is the fold, an integer between 1 and k
    # ... -> arguments to be passed to JGL
    # The function is parallelized for Windows
    omit <- lapply(all.folds, function(q) q[[i]]) # the indices of the observations in the ith fold in each class
    dat_train <- lapply(1:nclasses, function(cls) dat[[cls]][-omit[[cls]], , drop = FALSE]) # for each class, the data after omitting the ith fold = training set
    dat_test <- lapply(1:nclasses, function(cls) dat[[cls]][omit[[cls]], , drop = FALSE]) # for each class, the data after keeping only the ith fold = validation set
    Sigmas <- lapply(dat_test, cov) # for each class, the covariance matrices in the validation set
    
    if(ncores > 1)
    {
      cl <- makeCluster(ncores)
      clusterExport(cl, list("dat_train", "lambdas", "JGL"), envir=environment())
      Omegas <- parLapply(cl = cl, 1:nrow(lambdas), function(la) JGL(dat_train, lambda1 = lambdas[la, 1], lambda2 = lambdas[la, 2], return.whole.theta = TRUE, ...)$theta)
      stopCluster(cl)
      
    } else  {
      Omegas <- lapply(1:nrow(lambdas), function(la) JGL(dat_train, lambda1 = lambdas[la, 1], lambda2 = lambdas[la, 2], return.whole.theta = TRUE, ...)$theta) # for each pair of lambdas, the estimated partial correlation matrix in each class
    }
    
    fits <- sapply(Omegas, function(O) anpll(S = Sigmas, O = O))
    fits
  }
  
  # a matrix with the average negative predictive log likelihood function in each training sample
  anpllmat <- sapply(1:k, crossval,  all.folds = all.folds, nclasses = nclasses, ncores = ncores, lambdas = lambdas, ...)
  # average across folds
  anpllvec <- apply(anpllmat, 1, sum)
  
  # the pair of lambdas values with the minimum cost function in the validation set
  la <- lambdas[which.min(anpllvec), ]
  lambda1 <- unlist(la[1])
  lambda2 <- unlist(la[2])
  
  jgl <- JGL(dat, lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = TRUE, ...)
  
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





