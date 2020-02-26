dtp3 <-
  function(x, mu, par1, par2, FUN,param = "tp", log = FALSE ){
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
    ifelse( is.numeric(logPDF),ifelse( log, return(logPDF), return(exp(logPDF)) ), logPDF )
  }
