########################################################################################################
#
# DEFINITION OF CLASSES AND CORRESPONDING VALIDITY CHECKS
#
########################################################################################################


########################################################################################################
# CLASS normFit
########################################################################################################

setClass("normFit", representation("list"), prototype = prototype(elementType = "list"), contains="list")

valid_normFit <- function(object) {
  msg <- NULL
  if (!all(c('mu','Sigma','probs','probcluster','G') %in% names(object))) msg <- "normFit object missing either mu, Sigma, probs, probscluster or G"
  if (length(object$mu)!=object$G) msg <- "mu has the wrong length"
  if (length(object$Sigma)!=object$G) msg <- "Sigma has the wrong length"
  if (ncol(object$probs) != object$G) msg <- "probs has the wrong number of columns"
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("normFit", valid_normFit)



########################################################################################################
# CLASS skewtFit
########################################################################################################

setClass("skewtFit", representation("list"), prototype = prototype(elementType = "list"), contains="list")

valid_skewtFit <- function(object) {
  msg <- NULL
  if (!all(c('mu','Sigma','alpha','nu','probs','probcluster','G','ttype') %in% names(object))) msg <- "skewtFit object missing either mu, Sigma, alpha, nu, probs, probscluster or G"
  if (length(object$mu)!=object$G) msg <- "mu has the wrong length"
  if (length(object$Sigma)!=object$G) msg <- "Sigma has the wrong length"
  if (length(object$alpha)!=object$G) msg <- "alpha has the wrong length"
  if (length(object$nu)!=object$G) msg <- "nu has the wrong length"
  if (ncol(object$probs) != object$G) msg <- "probs has the wrong number of columns"
  if (!(object$ttype %in% c('independent','dependent'))) stop("ttype must be 'independent' or 'dependent'")
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("skewtFit", valid_skewtFit)
