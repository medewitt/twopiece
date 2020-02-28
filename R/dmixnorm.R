#' Density for mixture of Normals
#'
#' This function outputs the density for a mixture of normal distribution
#'
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken
#'     to be a quantile.
#' @param mu a list of means
#' @param Sigma a list covariance matrices
#' @param probs a numeric non-negative vector of length K, specifying the
#'     probability for the K classes
#' @param logscale return log density if \code{TRUE}, default \code{FALSE}
#'
#' @export
dmixnorm <- function(x, mu, Sigma, probs, logscale=FALSE) {
  K <- length(probs)

  if (length(mu) != K) {
    stop("mu must be a list with length = length(probs)")
  }

  if (length(Sigma) != K){
    stop("Sigma must be a list with length = length(probs)")
    }
  if (is.vector(x)){

    x <- matrix(x,nrow=1)
  }

  logd <- matrix(NA,nrow=nrow(x),ncol=K)

  for (i in 1:K) {
    logd[,i] <- mvtnorm::dmvnorm(x,mu[[i]], Sigma[[i]], log=TRUE) + log(probs[i])
    }
  ct <- apply(logd,1,max)

  ans <- exp(ct) * rowSums(exp(logd-ct))

  if (logscale) {
    ans <- log(ans)
  }

  return(ans)
}
