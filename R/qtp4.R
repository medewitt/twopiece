qtp4 <-
  function(p, mu, par1, par2, delta, FUN,param = "tp"){
    if(param == "tp")
    {
      ifelse( par1 > 0  & par2 > 0,
              Q  <-  ifelse( p < par1/(par1+par2) ,mu + par1*FUN(0.5*p*(par1+par2)/par1, delta),mu + par2*FUN(0.5*((par1+par2)*(1+p)-2*par1)/par2, delta)),
              Q  <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization tp' )
    }
    if(param == "eps")
    {
      sigma  = par1 ; gamma  = par2
      ifelse( sigma > 0 & abs(gamma) < 1,
              Q  <- ifelse( p < 0.5*(1+gamma),mu + sigma*(1+gamma)*FUN(p/(1+gamma), delta), mu + sigma*(1-gamma)*FUN((p-gamma)/(1-gamma), delta)) ,
              Q  <- 'invalid arguments: par1 is not positive or/and abs(par2) is not less than 1 in the parametrization eps' )
    }
    if(param == "isf")
    {
      sigma  = par1; gamma  = par2
      ifelse( sigma > 0 & gamma > 0,
              Q <-  ifelse( p < gamma^2/(1+gamma^2), mu + sigma*gamma*FUN(0.5*p*(1+gamma^2)/gamma^2, delta), mu + sigma*FUN(0.5*(p*(1+gamma^2)+1-gamma^2), delta)/gamma) ,
              Q <- 'invalid arguments: par1 or/and par2 is/are not positive in the parametrization isf' )

    }
    return(Q)
  }
