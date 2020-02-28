#Gibbs sampling for multivariate Normal mixture
# Input
# - x: observed data (observations in rows, variables in columns)
# - G: number of components
# - clusini: either vector with initial cluster allocation, 'kmedians' for K-medians (as in cclust from package flexclust) 'kmeans' for K-means, 'em' for EM algorithm based on mixture of normals (as in Mclust package, initialized by hierarchical clustering so it can be slow)
# - priorParam: list with named elements 'mu', 'Sigma', 'probs' containing prior parameters. See help(skewnormprior).
# - niter: number of Gibbs iterations
# - returnCluster: if set to TRUE the allocated cluster at each MCMC iteration is returned. This can be memory-consuming if nrow(x) is large.
# - relabel: 'none' to do no relabelling. 'ECR' for Papastamoulis-Iliopoulos 2010 to make simulated z similar to pivot z (taken from last iteration); 'RW' for Rodriguez-Walker (2014) relabelling (loss function aimed at preserving cluster means), 'PC' for identifiability constraint based on projection of location parameters on first principal component
# - verbose: set to TRUE to output iteration progress
# Output:
# - mu: list of length G, where mu[[i]] is a matrix with posterior draws (niter-burnin rows)
# - Sigma: list of length G, where Sigma[[i]] is a matrix with posterior draws (niter-burnin rows). Each row contains diag(S),S[upper.tri(S)]
# - probs: matrix with niter-burnin rows and G columns with posterior draws for the mixing probabilities
# - probcluster: matrix with nrow(x) rows and G columns with posterior probabilities that each observation belongs to each cluster (Rao-Blackwellized)
# - cluster: if returnCluster==TRUE, a matrix with niter-burnin rows and nrow(x) columns with latent cluster allocations at each MCMC iteration. If returnCluster==FALSE, NA is returned.

mixnormGibbs <- function(x, G, clusini='kmedians',
                         priorParam=normprior(ncol(x),G), niter,
                         burnin=round(niter/10), returnCluster=FALSE,
                         elabel='ECR', verbose=TRUE) {

  p <- ncol(x); n <- nrow(x)
  m <- priorParam[['mu']]$m; gprior <- priorParam[['mu']]$g
  Q <- priorParam[['Sigma']]$Q; q <- priorParam[['Sigma']]$q
  r <- priorParam[['probs']]['r']
  ## Initial cluster allocations ##
  if (G>1) {
    if (is.character(clusini)) {
      if (clusini=='kmedians') {
        z <- try(predict(flexclust::cclust(x, k=G, dist='manhattan', method='kmeans')))
      } else if (clusini=='kmeans') {
        z <- kmeans(x, centers=G)$cluster
      } else if (clusini=='em') {
        z <- try(mclust::Mclust(x, G=G, modelNames='VVV')$classification)
      } else stop("Invalid value for 'clusini'")
      if (class(z)=='try-error') z <- kmeans(x, centers=G)$cluster
    } else {
      if (length(clusini) != nrow(x)) stop("length(clusini) must be equal to nrow(x)")
      z <- as.integer(factor(clusini))
    }
  } else {
    z <- rep(1,nrow(x))
  }
  zcount <- rep(0,G); names(zcount) <- (1:G)
  tab <- table(z); zcount[names(tab)] <- tab
  ## Initialize parameters ##
  probscur <- as.vector(rdirichlet(1,r+zcount))
  mucur <- lapply(1:G, function(i) double(p))
  Scur <- lapply(1:G, function(i) matrix(NA,nrow=p,ncol=p))
  for (i in 1:G) {
    sel <- which(z==i)
    if (length(sel)>1) {
      mucur[[i]] <- apply(x[sel,,drop=FALSE],2,median); Scur[[i]] <- cov(x[sel,,drop=FALSE])
    } else {
      mucur[[i]] <- m; Scur[[i]] <- Q / (q+p+1)
    }
  }
  ##Gibbs sampling
  mu <- lapply(1:G, function(i) matrix(NA,nrow=niter-burnin, ncol=p))
  Sigma <- lapply(1:G, function(i) matrix(NA,nrow=niter-burnin, ncol=p*(p+1)/2))
  probs <- matrix(NA, nrow=niter-burnin, ncol=G)

  cluster <- matrix(NA,nrow=niter-burnin,ncol=nrow(x))
  colnames(cluster) <- paste('indiv',1:n,sep='')

  if (verbose) {
    niter10 <- round(niter/10)
    message("Running MCMC")
  }

  for (l in 1:niter) {
    #Sample latent cluster indicators z
    dx <- sapply(1:G,function(g) mvtnorm::dmvnorm(x,mucur[[g]],Scur[[g]],log=TRUE))
    dx <- t(t(dx)+log(probscur))
    dx <- exp(dx - rowMaxs(dx))
    #dx <- exp(dx - apply(dx,1,max))
    dx <- dx/rowSums(dx)
    z <- apply(dx,1,function(pp) match(TRUE,runif(1)<cumsum(pp)))

    tab <- table(z)
    zcount[1:G] <- 0
    zcount[names(tab)] <- tab

    #Sample component weights
    probscur <- as.vector(rdirichlet(1,alpha=r+zcount))

    #Sample mu, Sigma
    for (g in 1:G) {
      sel <- (z==g)
      if (sum(sel)>0) {
        newparam <- rmuSigmaNormGibbs(x[sel,,drop=FALSE],m=m,g=gprior,q=q,Q=Q)
        Scur[[g]] <- newparam$Sigma
        mucur[[g]] <- newparam$mu
      } else {
        #If no observations in cluster, sample from the prior
        Scur[[g]] <- rinvwishart(q,Q)
        mucur[[g]] <- as.vector(mvtnorm::rmvnorm(1,m,gprior * Scur[[g]]))
      }
    }
    #Save output
    if (l>burnin) {
      idx <- l-burnin
      for (g in 1:G) {
        mu[[g]][idx,] <- mucur[[g]]
        if (p>1) {
          Sigma[[g]][idx,] <- c(diag(Scur[[g]]),Scur[[g]][upper.tri(Scur[[g]])])
        } else {
          Sigma[[g]][idx,] <- Scur[[g]]
        }
      }
      probs[idx,] <- probscur
      cluster[idx,] <- z
    }
    if (verbose & ((l %% niter10)==0)) {
      cat('.')
      }
  }

  # Preparing Outputs
  names(mu) <- names(Sigma) <- colnames(probs) <- paste('cluster',1:G,sep='')

  ans <- list(mu=mu,
              Sigma=Sigma,
              probs=probs,
              probcluster=NA, cluster=NA,
              G=G)

  if (returnCluster) {
    ans$cluster <- cluster
    }

  fit <- new("normFit",ans)

  #Fix potential label switching issues

  if ((G>1) & (relabel!='none')) {
    if (verbose) {
      cat('\n');
      cat('Post-processing label-switching...')
      }

    fit <- fixLabelSwitch(fit,x=x,z=cluster,method=relabel)
    }
  #Compute cluster probabilities
  fit$probcluster <- clusterprobs(fit, x=x)
  rownames(fit$probcluster) <- paste('indiv',1:n,sep='')
  colnames(fit$probcluster) <- paste('cluster',1:G,sep='')

  if (verbose) {
    cat('\n')
  }

  return(fit)
}
