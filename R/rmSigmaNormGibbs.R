rmuSigmaNormGibbs <- function(x,m,g,q,Q) {
  #Gibbs sampler from the posterior of (mu,Sigma) given the latent cluster allocations
  #Input
  # - x: data (individuals in rows, variables in columns)
  # - m,g: prior on mu is N(m,g*Sigma)
  # - Q, q: prior on Sigma is InvWishart(Q,q), where q are the degrees of freedom
  #Output: list with 2 components
  # - mu: sampled value of mu
  # - Sigma: sampled value of Sigma
  n= nrow(x)
  if (n>1) { Spost= Q + cov(x) * (n-1) } else { Spost= Q }
  Sigma= rinvwishart(q+nrow(x), Spost)
  w= n/(n+1/g)
  mu= mvtnorm::rmvnorm(1, colMeans(x)*w + m*(1-w), Sigma/(n+1/g))
  return(list(mu=mu,Sigma=Sigma))
}
