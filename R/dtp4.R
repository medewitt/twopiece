dtp4 <-
  function(x, mu, par1, par2, delta, FUN,param = "tp", log = FALSE ){
    if(param == "tp")
    {
      ifelse( par1 > 0  & par2 > 0 & delta > 0,
              logPDF  <- log(2) + ifelse( x < mu ,FUN( (x-mu)/par1, delta, log=T), FUN( (x-mu)/par2, delta, log=T)) - log( par1+par2 ),
              logPDF  <- 'invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization tp' )
    }
    if(param == "eps")
    {
      sigma  = par1 ; gamma  = par2
      ifelse( sigma > 0 & abs(gamma) < 1 &  delta > 0,
              logPDF  <- ifelse( x < mu,FUN( (x-mu)/(sigma*(1+gamma)), delta, log=T ),FUN( (x-mu)/(sigma*(1-gamma)), delta, log=T ) ) - log( sigma ),
              logPDF  <- 'invalid arguments: par1 or/and delta is/are no positive or/and abs(par2) is no less that 1 in the parametrization eps' )
    }
    if(param == "isf")
    {
      sigma  = par1; gamma  = par2
      ifelse( sigma > 0 & gamma > 0  &  delta > 0,
              logPDF <- log(2) + ifelse( x < mu,FUN( (x-mu)/(sigma*gamma), delta, log=T ),FUN( (x-mu)/(sigma/gamma), delta, log=T) ) - log( sigma*(gamma+1/gamma) ),
              logPDF <- 'invalid arguments: par1 or/and par2 or/and delta is/are no positive in the parametrization isf' )

    }
    ifelse( is.numeric(logPDF),ifelse( log, return(logPDF), return(exp(logPDF)) ), logPDF )
  }
