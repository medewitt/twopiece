#' @aliases dtp3
#' @export
ptp3 <-
  function(x, mu, par1, par2, FUN, param = "tp", log.p = FALSE ){

    param = match.arg(param, choices = c("tp", "eps", "isf"))

    if(!is.logical(log.p)){
      stop("log.p must be a boolean")
    }

    if(param == "tp")
    {
      if(!(par1 > 0  & par2 > 0)){
        stop('invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp')
      }

      CDF  <-  ifelse( x < mu ,2*par1*FUN( (x-mu)/par1, log.p)/(par1+par2),
                       ( par1 + par2*(2*FUN( (x-mu)/par2, log.p)-1) )/( par1+par2))
    }
    if(param == "eps")
    {
      sigma  = par1
      gamma  = par2

      if(!(sigma > 0 & abs(gamma) < 1)){
        stop('invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps')
      }

      CDF  <- ifelse( x < mu,
                      (1+gamma)*FUN( (x-mu)/(sigma*(1+gamma)), log.p ),
                      gamma + (1-gamma)*FUN( (x-mu)/(sigma*(1-gamma)), log.p ) )
    }
    if(param == "isf")
    {
      if(!(sigma > 0 & gamma > 0)){
        stop('invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf')
      }
      sigma  = par1; gamma  = par2

      CDF <-  ifelse( x < mu,
                      2*gamma^2*FUN( (x-mu)/(sigma*gamma), log.p=F )/(1+gamma^2),
                      ( gamma^2-1 + 2*FUN( (x-mu)/(sigma/gamma), log.p)  )/(1+gamma^2))

    }
    ifelse( log.p,
            return(log(CDF)),
            return(CDF) )
  }
