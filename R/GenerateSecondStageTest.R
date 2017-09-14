
#' Title
#'
#' @param delta
#' @param sigma
#' @param M
#' @param second.stage.mat
#'
#' @return
#' @export
#'
#' @examples
GenerateSecondStageTest <- function(delta,sigma,M,second.stage.mat){
  delta.no.inter <- delta[(M+1):length(delta)]
  ind.delta <- 0
  ind.covar <- 0
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

  all.covar.names <- colnames(second.stage.mat)
  all.second.stage.names <- row.names(second.stage.mat)

  covar.names <- NULL
  second.stage.effect.names <- NULL
  for(i in 1:nrow(second.stage.mat)){
    for(j in 1:ncol(second.stage.mat)){
      if(is.na(second.stage.mat[i,j])==F){
        covar.names = c(covar.names,all.covar.names[j])
        second.stage.effect.names =
          c(second.stage.effect.names,all.second.stage.names[i])
      }
    }
  }




  result <- data.frame(CovarName = covar.names,
                       SecondStageEffect = second.stage.effect.names,
                       odds,
                       odds.low,
                       odds.high,
                       p.individual.heter)

  colnames(result) <- c("CovarName","SecondStageEffect",
                        "OddsRatio",
                        "OddsRatio(95%CI low)",
                        "OddsRatio(95%CI high)",
                        "Pvalue")
  return(result)

}
