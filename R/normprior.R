#Create list containing prior parameters for multivariate Normal mixture. The following prior is used:
#  (mu,Sigma) ~ N(mu; m,g*Sigma) * IWishart(Sigma; Q, q)
#  probs ~ Symmetric Dirichlet(r), where r is a scalar
# Input
# - p: dimension of the observed data
# - G: number of components
# - m, g: parameters for prior on mu given Sigma
# - Q, q: parameters for prior on Sigma
# - r: parameter for symmetric Dirichlet prior on mixing probabilities

normprior <- function(p,G,m=rep(0,p),g=1,Q=diag(p),q=p+1,r=1/G) {
  if (length(m) != p) stop("m has the wrong length")
  if (length(g)>1) stop("g must have length 1")
  if (!is.matrix(Q)) stop("Q must be a matrix")
  if ((nrow(Q) != ncol(Q)) | (nrow(Q) != p)) stop("Q must be a square matrix with p rows")
  if (any(eigen(Q,symmetric=TRUE)$values <= 0)) stop("Q is not positive-definite!")
  if (q < p) stop("q must be >= p")
  if (length(r)>1) stop("r must have length 1")
  ans <- vector("list",5)
  names(ans) <- c('mu','Sigma','probs')
  ans[['mu']] <- list(m=m,g=g)
  ans[['Sigma']] <- list(Q=Q,q=q)
  ans[['probs']] <- c(r=r)
  return(ans)
}
