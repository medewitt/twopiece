#' @title
#' The 3-Parameter Two Piece Distribution
#' @description
#' Density, distribution function, quantile function and random
#' generation for the 3-parameter two piece distribution with 3
#' parameterizations: two-piece (tp), epsilon-skew (eps), and inverse scale factors (isf).
#'
#' @param x vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu  location parameter, \eqn{\mu}
#' @param par1 scale parameter 1, \eqn{\sigma_1}
#' @param par2 scale parameter 2, \eqn{\sigma_2}
#' @param FUN  a symmetric density \code{f}
#' @param param parameterizations used with default of \code{"tp"}
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#'
#' @usage
#' dtp3(x, mu, par1, par2, FUN, param = "tp", log = FALSE)
#' ptp3(x, mu, par1, par2, FUN, param = "tp", log.p = FALSE)
#' qtp3(p, mu, par1, par2, FUN, param = "tp")
#' rtp3(n, mu, par1, par2, FUN, param = "tp")
#'
#' @details
#' The 3-parameter two piece distribution with parameters \eqn{\mu}
#' \eqn{\sigma_1}and \eqn{\sigma_2}
#' has the following density:
#'
#' \deqn{s(x) = \frac{2}{\sigma_1 + \sigma_2}f( ( x -\mu )/\sigma_1 ) \quad for \quad x<\mu }
#' and
#' \deqn{s(x) = \frac{2}{\sigma_1 + \sigma_2}f( ( x -\mu )/\sigma_2 ) \quad for \quad x\geq \mu }
#'
#' where \code{f(x)} is a symmetric density about zero.
#'
#' If param is not specified, it assumes the default value of "tp".
#' Information about the "eps" and "isf" parameterizations can be found in the References.
#' \code{dtp3} gives the density, \code{ptp3} gives the distribution function,
#' \code{qtp3} gives the quantile function and \code{rtp3} generates random deviates.
#'
#'
#' @author
#' F. J. Rubio \email{fxrubio@gmail.com} and A. M. López \email{amonteslop@gmail.com}
#'
#' @references
#' Arellano-Valle, R. B Gómez, H. W. and Quintana, F. A. (2005). Statistical inference for general class of
#'asymmetric distributions. \emph{Journal of Statistical Planning and Inference}, \bold{128}: 427-443.
#'
#' Fernández, C. and Steel, M. F. J. (1998). On Bayesian modeling of fat tails and skewness.
#' \emph{Journal of the American Statistical Asociation}, \bold{93}, 359-371.
#'
#' Mudholkar, G. S. and Hutson, A. D. (2000). The epsilon-skew-normal distribution for analyzing
#' near-normal data. \emph{Journal of Statistical Planning and Inference}, \bold{83}: 291-309.
#'
#' Rubio, F. J. and Steel, M. F. J. (2014). Inference in Two-Piece Location-Scale models with
#' Jeffreys Priors, with discussion. \emph{Bayesian Analysis}, \bold{9}: 1-22.
#'
#' @seealso
#' \link{dnorm} for the normal distribution and \link{dt} for the Student t distribution
#' \code{\link{dtp4}} for the 4-parameter two piece distribution.
#'
#'
#' @examples
#'
#' ## 3-parameter two piece normal density with parameterization "tp"
#' tempf = function(x) dtp3(x,0,3,1,dnorm,param="tp")
#' curve(tempf,-10,5)
#'
#' ## 3-parameter two piece normal distribution with parameterization "tp"
#' tempf = function(x) ptp3(x,0,1,3,pnorm,param="tp")
#' curve(tempf,-10,10)
#'
#' ## random number generation for 3-parameter two piece normal distribution
#' ## with parameterization "tp"
#' sim <- rtp3(1000,0,1,1,rnorm)
#' hist(sim,probability=TRUE)
#'
#' ## quantile function for the 3-parameter two piece normal distribution
#' ## with parameterization "tp"
#' qtp3(0.5, 0, 1, 1, qnorm, param = "tp")
#' @name dtp3
#' @aliases dtp3
#' @aliases qtp3
#' @aliases rtp3
#' @aliases ptp3
#' @export
#'

dtp3 <-
  function(x, mu, par1, par2, FUN, param = "tp", log = FALSE ){

    param = match.arg(param, choices = c("tp", "eps", "isf"))

    if(param == "tp")
    {
      ifelse( par1 > 0  & par2 > 0,
              logPDF  <- log(2) + ifelse( x < mu ,FUN( (x-mu)/par1, log=T), FUN( (x-mu)/par2, log=T)) - log( par1+par2 ),
              logPDF  <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp' )
    }
    if(param == "eps")
    {
      sigma  = par1 ; gamma  = par2
      ifelse( sigma > 0 & abs(gamma) < 1,
              logPDF  <- ifelse( x < mu,FUN( (x-mu)/(sigma*(1+gamma)), log=T ),FUN( (x-mu)/(sigma*(1-gamma)), log=T ) ) - log( sigma ),
              logPDF  <- 'invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps' )
    }
    if(param == "isf")
    {
      sigma  = par1; gamma  = par2
      ifelse( sigma > 0 & gamma > 0,
              logPDF <- log(2) + ifelse( x < mu,FUN( (x-mu)/(sigma*gamma), log=T ),FUN( (x-mu)/(sigma/gamma), log=T) ) - log( sigma*(gamma+1/gamma) ),
              logPDF <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf' )

    }
    ifelse( is.numeric(logPDF),
            ifelse( log, return(logPDF),
                    return(exp(logPDF)) ),
            logPDF )
  }
