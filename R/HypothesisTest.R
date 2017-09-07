#' Title
#'
#' @param logodds
#' @param infor
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
  p.value.GTA <- format(p.value.GTA,scientific = T,digits = 3)
  return(p.value.GTA)

}


#' Title
#'
#' @param logodds
#' @param infor
#'
#' @return
#' @export
#'
#' @examples
GlobalTestForHeter <- function(logodds,sigma){

  sigma <- as.matrix(sigma)
  df <- length(logodds)
  sigma.casecase <- sigma[2:df,2:df]
  logodds.casecase <- logodds[2:df]
  GTH.stat <- t(logodds.casecase)%*%solve(sigma.casecase)%*%logodds.casecase
  p.value.GTH <- pchisq(as.numeric(GTH.stat),df=(df-1),lower.tail = F)
  p.value.GTH <- format(p.value.GTH,scientific = T,digits = 3)
  return(p.value.GTH)
}

#' Title
#'
#' @param logodds
#' @param infor
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
  p.value.IHT <- format(p.value.IHT,scientific = T, digits=3)
  return(p.value.IHT)
}
