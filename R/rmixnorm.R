#' Multivariate draws from Normal mixture
#'
#' This function generates a specified number of draws from
#' a mixture of normal distributions
#'
#' @param n the number of draws
#' @param mu a list of means
#' @param Sigma a list of convariance matrices for the distributions
#' @param probs a numeric non-negative vector of length K, specifying
#'     the probability for the K classes
#'
#' @export
#'
rmixnorm <- function(n,mu,Sigma,probs) {
  #Checks
  K <- length(probs)

  if (!is.list(mu)) {

    mu <- lapply(1:K,function(z) mu)

    warning("Formatted mu as list")

    } else {
      if (length(mu) != K) {
        stop("mu has the wrong length")
      }
    }
  if (!is.list(Sigma)) {

    Sigma <- lapply(1:K, function(z) Sigma);

    warning("Formatted Sigma as list")

    } else {
      if (length(Sigma) != K) {
        stop("Sigma has the wrong length")
      }
    }
  p <- length(mu[[1]])

  #Obtain draws
  ngroup <- as.vector(rmultinom(1, size=n, prob=probs))

  groupidx <- c(0, cumsum(ngroup))

  x <- matrix(NA, nrow=n, ncol=p)

  for (i in 1:K) {
    x[(groupidx[i]+1):groupidx[i+1],] <- mvtnorm::rmvnorm(ngroup[i],mu[[i]],Sigma[[i]])
  }

  idx <- sample(1:nrow(x), size=nrow(x), replace=FALSE)

  x <- x[idx,]

  cluster <- rep(1:K,ngroup)[idx]

  return(list(x=x,
              cluster=cluster))
}
