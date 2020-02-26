################################################################################################################
## METHODS FOR OBJECTS OF CLASS normFit
################################################################################################################

setMethod("show", signature(object="normFit"), function(object) {
  cat("Normal mixture with ",object$G," components\n")
  cat("Use coef() to obtain posterior means. clusterprobs() for cluster probabilities\n")
}
)

setMethod("coef", signature(object="normFit"), function(object, ...) {
  mu <- lapply(object$mu,colMeans)
  p <- length(mu[[1]])
  diag <- 1:p; if (p>1) { nondiag <- (p+1):ncol(object$Sigma[[1]]) } else { nondiag= integer(0) }
  Sigma <- lapply(object$Sigma, function(z) { vec2matrix(colMeans(z),diag=diag,nondiag=nondiag) })
  probs <- colMeans(object$probs)
  return(list(mu=mu,Sigma=Sigma,probs=probs))
}
)


setGeneric("coefMedian", function(object,...) standardGeneric("coefMedian"))

setMethod("coefMedian", signature(object="normFit"), function(object, ...) {
  colMedians <- function(z,...) apply(z,2,'median',...)
  mu <- lapply(object$mu,colMedians)
  p <- length(mu[[1]])
  diag <- 1:p; if (p>1) { nondiag <- (p+1):ncol(object$Sigma[[1]]) } else { nondiag= integer(0) }
  Sigma <- lapply(object$Sigma, function(z) { vec2matrix(colMedians(z),diag=diag,nondiag=nondiag) })
  probs <- colMedians(object$probs)
  return(list(mu=mu,Sigma=Sigma,probs=probs))
}
)


#Compute cluster probabilities under Normal mixture, averaging over posterior draws contained in x
# - fit: object of type normFit
# - x: data points at which to evaluate the probabilities
# - iter: iterations of posterior draws over which to average. Defaults to using 1,000 equally spaced iterations.
setMethod("clusterprobs", signature(fit='normFit'), function(fit, x, iter) {
  if (missing(iter)) {
    if (nrow(fit$probs)>1000) {
      iter <- as.integer(seq(1,nrow(fit$probs),length=min(nrow(fit$probs),1000)))
    } else {
      iter <- 1:nrow(fit$probs)
    }
  }
  if (!is.integer(iter)) stop("iter must contain integer iteration numbers")
  if ((min(iter)<1) | max(iter)>nrow(fit$mu[[1]])) stop("Specified iter values are out of valid range")
  if (missing(x)) stop("x values must be provided")
  G <- fit$G; p <- ncol(fit$mu[[1]])
  if (p>1) {
    diag <- 1:p
    nondiag <- (p+1):ncol(fit$Sigma[[1]])
    S <- lapply(1:G, function(i) { ans <- apply(fit$Sigma[[i]][iter,,], 1, vec2matrix, diag=diag, nondiag=nondiag); array(as.vector(ans),dim=c(p,p,ncol(ans))) } )
    proboneiter <- function(i) {
      dx <- sapply(1:G,function(g) dmvnorm(x,fit$mu[[g]][iter[i],],S[[g]][,,i],log=TRUE))
      dx <- t(t(dx)+log(fit$probs[i,]))
      dx <- exp(dx - rowMaxs(dx))
      dx/rowSums(dx)
    }
    ans <- lapply(1:length(iter),function(i) proboneiter(i))
  } else {
    proboneiteruniv <- function(i) {
      dx <- sapply(1:G,function(g) dmvnorm(x,fit$mu[[g]][iter[i],],fit$Sigma[[g]][iter[i],,drop=FALSE],log=TRUE))
      dx <- t(t(dx)+log(fit$probs[i,]))
      dx <- exp(dx - rowMaxs(dx))
      dx/rowSums(dx)
    }
    ans <- lapply(1:length(iter),function(i) proboneiteruniv(i))
  }
  ans <- Reduce('+',ans)/length(ans)
  return(ans)
}
)

loglmeanNorm <- function(fit,x) {
  #log-likelihood evaluated at posterior mean
  parest <- coef(fit)
  sum(dmixnorm(x,mu=parest$mu,Sigma=parest$Sigma,probs=parest$probs,logscale=TRUE))
}

loglNorm <- function(fit,x) {
  #log-likelihood at each MCMC iteration
  if (class(fit) != 'normFit') stop("fit must be of class normFit")
  G <- fit$G; p <- ncol(fit$mu[[1]])
  diag <- 1:p; if (p>1) { nondiag <- (p+1):ncol(fit$Sigma[[1]]) } else { nondiag= integer(0) }
  S <- lapply(1:G, function(i) { ans <- apply(fit$Sigma[[i]], 1, vec2matrix, diag=diag, nondiag=nondiag); array(as.vector(ans),dim=c(p,p,ncol(ans))) } )
  lhoodoneiter <- function(i) {
    dx <- sapply(1:G,function(g) dmvnorm(x,fit$mu[[g]][i,],S[[g]][,,i],log=TRUE))
    dx <- t(t(dx)+log(fit$probs[i,]))
    sum(log(rowSums(exp(dx))))
  }
  ans <- sapply(1:nrow(fit$mu[[1]]),function(i) lhoodoneiter(i))
  ans
}


setMethod("fixLabelSwitch", signature(fit="normFit"), function(fit,x,z,method='ECR') {
  #Permute component labels to avoid label switching issues. Components relabelled to have increasing projection on first PC
  # Input
  # - fit: object of type normFit
  # - x: observed data used to obtain fit
  # - z: latent cluster allocations at each MCMC iteration
  # - method: 'ECR' for Papastamoulis-Iliopoulos 2010 to make simulated z similar to pivot z (taken from last iteration); 'RW' for Rodriguez-Walker (2014) relabelling (loss function aimed at preserving cluster means), 'PC' for identifiability constraint based on projection of location parameters on first principal component
  # Output: object of type normFit with relabelled components
  if (method=='ECR') {
    r <- ecr(z[nrow(z),],z=z,K=fit$G)[[1]]
  } else if (method=='RW') {
    r <- dataBased(x,K=fit$G,z)[[1]]
  } else if (method=='PC') {
    e <- eigen(cov(x))$vectors[,1,drop=FALSE]
    proj <- do.call(cbind,lapply(fit$mu, "%*%", e))
    r <- t(apply(proj,1,rank,ties.method='first'))
  } else { stop("Invalid value for 'method' in fixLabelSwitch") }
  fitnew <- fit
  for (g in 1:fit$G) {
    glabel <- apply(r==g, 1, function(z) match(TRUE,z)) #component in fit corresponding to g^th reordered component in fitnew
    for (gg in 1:fit$G) {
      sel <- glabel==gg
      fitnew$mu[[g]][sel,] <- fit$mu[[gg]][sel,]
      fitnew$Sigma[[g]][sel,] <- fit$Sigma[[gg]][sel,]
      fitnew$probs[sel,g] <- fit$probs[sel,gg]
      if (!is.null(dim(fit$cluster))) fitnew$cluster[fit$cluster==gg] <- g
    }
  }
  return(fitnew)
}
)



################################################################################################################
## ROUTINES TO EVALUATE DENSITY, RANDOM NUMBER GENERATION
################################################################################################################


#Density for mixture of Normals (arguments as in dmvnorm)
dmixnorm <- function(x,mu,Sigma,probs,logscale=FALSE) {
  K <- length(probs)
  if (length(mu) != K) stop("mu must be a list with length = length(probs)")
  if (length(Sigma) != K) stop("Sigma must be a list with length = length(probs)")
  if (is.vector(x)) x <- matrix(x,nrow=1)
  logd <- matrix(NA,nrow=nrow(x),ncol=K)
  for (i in 1:K) { logd[,i] <- dmvnorm(x,mu[[i]],Sigma[[i]],log=TRUE) + log(probs[i]) }
  ct <- apply(logd,1,max)
  ans <- exp(ct) * rowSums(exp(logd-ct))
  if (logscale) ans <- log(ans)
  return(ans)
}



#Multivariate draws from Normal mixture
rmixnorm <- function(n,mu,Sigma,probs) {
  #Checks
  K <- length(probs)
  if (!is.list(mu)) { mu <- lapply(1:K,function(z) mu); warning("Formatted mu as list") } else { if (length(mu) != K) stop("mu has the wrong length") }
  if (!is.list(Sigma)) { Sigma <- lapply(1:K,function(z) Sigma); warning("Formatted Sigma as list") } else { if (length(Sigma) != K) stop("Sigma has the wrong length") }
  p <- length(mu[[1]])
  #Obtain draws
  ngroup <- as.vector(rmultinom(1, size=n, prob=probs))
  groupidx <- c(0,cumsum(ngroup))
  x <- matrix(NA,nrow=n,ncol=p)
  for (i in 1:K) { x[(groupidx[i]+1):groupidx[i+1],] <- rmvnorm(ngroup[i],mu[[i]],Sigma[[i]]) }
  idx <- sample(1:nrow(x), size=nrow(x), replace=FALSE)
  x <- x[idx,]
  cluster <- rep(1:K,ngroup)[idx]
  return(list(x=x,cluster=cluster))
}



################################################################################################################
## GIBBS SAMPLING
################################################################################################################

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
mixnormGibbs <- function(x, G, clusini='kmedians', priorParam=normprior(ncol(x),G), niter, burnin=round(niter/10), returnCluster=FALSE, relabel='ECR', verbose=TRUE) {
  p <- ncol(x); n <- nrow(x)
  m <- priorParam[['mu']]$m; gprior <- priorParam[['mu']]$g
  Q <- priorParam[['Sigma']]$Q; q <- priorParam[['Sigma']]$q
  r <- priorParam[['probs']]['r']
  ## Initial cluster allocations ##
  if (G>1) {
    if (is.character(clusini)) {
      if (clusini=='kmedians') {
        z <- try(predict(cclust(x, k=G, dist='manhattan', method='kmeans')))
      } else if (clusini=='kmeans') {
        z <- kmeans(x, centers=G)$cluster
      } else if (clusini=='em') {
        z <- try(Mclust(x, G=G, modelNames='VVV')$classification)
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
  mu <- lapply(1:G,function(i) matrix(NA,nrow=niter-burnin,ncol=p))
  Sigma <- lapply(1:G,function(i) matrix(NA,nrow=niter-burnin,ncol=p*(p+1)/2))
  probs <- matrix(NA,nrow=niter-burnin,ncol=G)
  cluster <- matrix(NA,nrow=niter-burnin,ncol=nrow(x)); colnames(cluster) <- paste('indiv',1:n,sep='')
  if (verbose) { niter10 <- round(niter/10); cat("Running MCMC") }
  for (l in 1:niter) {
    #Sample latent cluster indicators z
    dx <- sapply(1:G,function(g) dmvnorm(x,mucur[[g]],Scur[[g]],log=TRUE))
    dx <- t(t(dx)+log(probscur))
    dx <- exp(dx - rowMaxs(dx))
    #dx <- exp(dx - apply(dx,1,max))
    dx <- dx/rowSums(dx)
    z <- apply(dx,1,function(pp) match(TRUE,runif(1)<cumsum(pp)))
    tab <- table(z); zcount[1:G] <- 0; zcount[names(tab)] <- tab
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
        mucur[[g]] <- as.vector(rmvnorm(1,m,gprior * Scur[[g]]))
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
    if (verbose & ((l %% niter10)==0)) cat('.')
  }
  names(mu) <- names(Sigma) <- colnames(probs) <- paste('cluster',1:G,sep='')
  ans <- list(mu=mu,Sigma=Sigma,probs=probs,probcluster=NA,cluster=NA,G=G)
  if (returnCluster) ans$cluster <- cluster
  fit <- new("normFit",ans)
  #Fix potential label switching issues
  if ((G>1) & (relabel!='none')) {
    if (verbose) { cat('\n'); cat('Post-processing label-switching...') }
    fit <- fixLabelSwitch(fit,x=x,z=cluster,method=relabel)
  }
  #Compute cluster probabilities
  fit$probcluster <- clusterprobs(fit, x=x)
  rownames(fit$probcluster) <- paste('indiv',1:n,sep='')
  colnames(fit$probcluster) <- paste('cluster',1:G,sep='')
  if (verbose) cat('\n')
  return(fit)
}



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
  mu= rmvnorm(1, colMeans(x)*w + m*(1-w), Sigma/(n+1/g))
  return(list(mu=mu,Sigma=Sigma))
}
