#' @aliases dtp3
#' @export
#'
rtp3 <-
  function(n, mu, par1, par2, FUN, param = "tp"){

    param = match.arg(param, choices = c("tp", "eps", "isf"))

    if(param == "tp")
    {
      ifelse( par1 > 0  & par2 > 0,
              sample  <- ifelse(runif(n)<par1/(par1+par2), mu - par1*abs(FUN(n)), mu + par2*abs(FUN(n))),
              sample  <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp' )
    }
    if(param == "eps")
    {
      sigma  = par1 ; gamma  = par2
      ifelse( sigma > 0 & abs(gamma) < 1,
              sample  <- ifelse(runif(n)<0.5*(1+gamma), mu - sigma*(1+gamma)*abs(FUN(n)), mu + sigma*(1-gamma)*abs(FUN(n))),
              sample  <- 'invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps' )
    }
    if(param == "isf")
    {
      sigma  = par1; gamma  = par2
      ifelse( sigma > 0 & gamma > 0,
              sample  <- ifelse(runif(n)<gamma^2/(1+gamma^2), mu - sigma*gamma*abs(FUN(n)), mu + sigma*abs(FUN(n))/gamma),
              sample  <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf' )

    }
    return(sample)
  }
