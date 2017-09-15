
#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
GlobalTestForAssoc <- function(logodds,sigma){
  sigma <- as.matrix(sigma)
  df <- length(logodds)
  GTA.stat <- t(logodds)%*%solve(sigma)%*%logodds
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)

  return(p.value.GTA)

}



#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
GlobalTestForHeter <- function(logodds,sigma){

  if(length(logodds)==1){
    return(NA)
  }else{
    sigma <- as.matrix(sigma)
    df <- length(logodds)
    sigma.casecase <- sigma[2:df,2:df]
    logodds.casecase <- logodds[2:df]
    GTH.stat <- t(logodds.casecase)%*%solve(sigma.casecase)%*%logodds.casecase
    p.value.GTH <- pchisq(as.numeric(GTH.stat),df=(df-1),lower.tail = F)
    places <- 3
    power.number <- floor(-log10(p.value.GTH))+places
    p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)


    return(p.value.GTH)
  }


}


#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
IndividualHeterTest <- function(logodds,sigma){

  sigma <- as.matrix(sigma)
  var.logodds <- diag(sigma)
  df <- length(logodds)
  z <- logodds/sqrt(var.logodds)
  p.value.IHT <- PvalueFunction(z)

  return(p.value.IHT)

}
