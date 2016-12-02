

#   Patrick Danaher's code (email 25/04/2016)
#   ROCmat[r,22] = n[1]*(sum(temp$theta[[1]]*S1) - log(det(temp$theta[[1]]))) +
#     n[2]*(sum(temp$theta[[2]]*S2) - log(det(temp$theta[[2]]))) +
#     n[3]*(sum(temp$theta[[3]]*S3) - log(det(temp$theta[[3]]))) +
#     2*(ROCmat[r,3] + ROCmat[r,4])		#AIC    
#   ROCmat[r,23] = sum(temp$theta[[1]]*S1) - log(det(temp$theta[[1]])) +
#     sum(temp$theta[[2]]*S2) - log(det(temp$theta[[2]])) +
#     sum(temp$theta[[3]]*S3) - log(det(temp$theta[[3]])) +
#     log(n[1])*(ROCmat[r,3] + ROCmat[r,4])	#BIC      (this is what Guo used)


AIC_jgl <- function(jgl, n, S, ...)
{
  # jgl = output of JGL
  # n = vector of sample sizes corresponding to each element of jgl (same order)
  # S = list of covariance matrices corresponding to each element of jgl (same order)
  
  # Formula 6.21 in Danaher et al. (2014) doi:10.1111/rssb.12033
  #   # notice that sum(jgl$theta[[i]] * S[[i]]) is the same as sum(diag(jgl$theta[[i]] %*% S[[i]])) but requires much less computational time
  aics <- sapply(1:length(jgl$theta),
                 function(k) n[k] * sum(S[[k]]*jgl$theta[[k]]) - n[k] * log(det(jgl$theta[[k]])) +
                   2*sum(jgl$theta[[k]][lower.tri(jgl$theta[[k]])] != 0))
  sum(aics)
  
  # altenrative code  
  #   logli <- sapply(1:length(jgl$theta),
  #           function(k)  n[k] * .5 * (log(det(jgl$theta[[k]])) - sum(S[[k]]*jgl$theta[[k]])))
  #   
  #   npar <- sapply(1:length(jgl$theta),
  #          function(k) sum(jgl$theta[[k]][lower.tri(jgl$theta[[k]])] != 0))
  # c("AIC" = sum(2*npar - 2*logli), "logli" = sum(logli), "npar" = sum(npar))
}


AIC_jgl2 <- function(jgl, n, S, dec = 5, ...)
{
  # jgl = output of JGL
  # n = vector of sample sizes corresponding to each element of jgl (same order)
  # S = list of covariance matrices corresponding to each element of jgl (same order)
  # dec = edges are rounded to the dec decimal place when considering them different for computing the AIC
  
  # Formula 6.21 in Danaher et al. (2014) doi:10.1111/rssb.12033
  # notice that sum(jgl$theta[[i]] * S[[i]]) is the same as sum(diag(jgl$theta[[i]] %*% S[[i]])) but requires much less computational time
  # in AIC_jgl2 the parameters that are equal across classes under a certain tol value are counted only once and not once by class
  L <- jgl$theta
  arr <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  arr <- round(arr, dec)
  Nuni <- apply(arr, 1:2, function(x) {un <- unique(x); length(un[un!=0])})
  
  aics <- sapply(1:length(jgl$theta),
                 function(k) n[k] * sum(S[[k]]*jgl$theta[[k]]) - n[k] * log(det(jgl$theta[[k]]))) 
  sum(aics) + 2*sum(Nuni[lower.tri(Nuni)])
}

logGaus <- function(S,K,n)
{
  KS = K %*% S
  tr = function(A) sum(diag(A))
  return(n/2 * (log(det(K)) - tr(KS))  )
}

# Computes the EBIC:
# EBIC <- function(S,K,n,gamma = 0.5,E,countDiagonal=FALSE)
# {
#   #   browser()
#   L <- logGaus(S, K, n)
#   if (missing(E)){
#     E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
#   }
#   p <- nrow(K)
#   
#   # return EBIC:
#   -2 * L + E * log(n) + 4 * E * gamma * log(p)  
# }

BIC_jgl <- function(jgl, n, S, gamma = 0, dec = 5)
{
  # This is the formula in Guo et al., 2011
  # jgl = output of JGL
  # n = vector of sample sizes corresponding to each element of jgl (same order)
  # S = list of covariance matrices corresponding to each element of jgl (same order)

  # Number of unique pars:
  L <- jgl$theta
  arr <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  arr <- round(arr, dec)
  Nuni <- apply(arr, 1:2, function(x) {
    un <- unique(x)
    length(un[un!=0])
    })
                
  # Number of pars:
  E <- sum(Nuni[upper.tri(Nuni)])
  
  # Number of nodes:
  p <- ncol(S[[1]])
  
  # likelihoods
  Lik = numeric(length(S))
  EBIC = numeric(length(S))
  for (i in seq_along(S)){
    Lik[i] = logGaus(S[[i]], jgl$theta[[i]], n[i])
  }
  
  # Sum likelihoods:
  L <- sum(Lik)
  
  # Get EBIC:
  -2 * L + E * log(sum(n)) + 4 * E * gamma * log(p)  
}
