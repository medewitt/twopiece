#' @export
ptp4 <-
  function(x, mu, par1, par2, delta, FUN,param = "tp", log.p = FALSE ){

    param = match.arg(param, choices = c("tp", "eps", "isf"))

    if(!is.logical(log.p)){
      stop("log.p must be a boolean")
    }

    if(param == "tp")
    {
      ifelse( par1 > 0  & par2 > 0,
              CDF  <-  ifelse( x < mu ,2*par1*FUN( (x-mu)/par1, delta, log.p = F)/(par1+par2),( par1 + par2*(2*FUN( (x-mu)/par2, delta, log.p = F)-1) )/( par1+par2)),
              CDF  <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp' )
    }
    if(param == "eps")
    {
      sigma  = par1
      gamma  = par2
      ifelse( sigma > 0 & abs(gamma) < 1,
              CDF  <- ifelse( x < mu, (1+gamma)*FUN( (x-mu)/(sigma*(1+gamma)), delta, log.p  = F),gamma + (1-gamma)*FUN( (x-mu)/(sigma*(1-gamma)), delta, log.p  = F) ) ,
              CDF  <- 'invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps' )
    }
    if(param == "isf")
    {
      sigma  = par1
      gamma  = par2
      ifelse( sigma > 0 & gamma > 0,
              CDF <-  ifelse( x < mu,2*gamma^2*FUN( (x-mu)/(sigma*gamma), delta, log.p = F)/(1+gamma^2),  ( gamma^2-1 + 2*FUN( (x-mu)/(sigma/gamma), delta, log.p = F)  )/(1+gamma^2)) ,
              CDF <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf' )

    }
    ifelse( is.numeric(CDF),ifelse( log.p, return(log(CDF)), return(CDF) ), CDF )
  }
