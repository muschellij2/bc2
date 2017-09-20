
#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples

DisplayTestResult = function(logodds,sigma){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  places <- 2
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.global.assoc <- GlobalTestForAssoc(logodds,sigma)
  p.global.heter <- GlobalTestForHeter(logodds,sigma)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = data.frame(matrix(0,1,2*length(odds)+2))
  for(i in 1:length(logodds)){
    result[1,2*i-1] <- paste0(odds[i],"(",odds.low[i],"-",
                             odds.high[i],")")
    result[1,2*i] <- p.individual.heter[i]
  }
  result[,2*length(odds)+1] <- p.global.assoc
  result[,2*length(odds)+2] <- p.global.heter
  return(result)
}



DisplayIndTestResult = function(logodds,sigma){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  places <- 2
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = NULL
  for(i in 1:length(logodds)){
    temp <- c(odds,odds.low,odds.high,p.individual.heter[i])
    result= rbind(result,temp)
  }
  return(result)
}
