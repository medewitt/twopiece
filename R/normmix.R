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







################################################################################################################
## GIBBS SAMPLING
################################################################################################################










