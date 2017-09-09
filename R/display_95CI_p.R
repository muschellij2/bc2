
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
  places <- 3
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.global.assoc <- GlobalTestForAssoc(logodds,sigma)
  p.global.heter <- GlobalTestForHeter(logodds,sigma)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = NULL
  for(i in 1:length(logodds)){
    result= c(result,paste0(odds[i],"(",odds.low[i],"-",
                            odds.high[i],")"),
              p.individual.heter[i])
  }
  result <- c(result,p.global.assoc, p.global.heter)
  return(result)
}



DisplayIndTestResult = function(logodds,sigma){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = NULL
  for(i in 1:length(logodds)){
    result= c(result,paste0(logodds[i],"(",logodds.low[i],"-",
                            logodds.high[i],")"),
              p.individual.heter[i])
  }
  return(result)
}
