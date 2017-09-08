#' Title
#'
#' @param baselineonly.number
#' @param main.effect.number
#' @param pairwise.interaction.number
#' @param saturated.number
#' @param baselineonly.second.cat
#' @param main.effect.second.cat
#' @param pairwise.interaction.second.cat
#' @param saturated.second.cat
#' @param M
#' @param total.covar.number
#' @param full.second.stage.names
#' @param covar.names
#'
#' @return
#' @export
#'
#' @examples
GenerateSecondStageMat <- function(baselineonly.number,
                                   main.effect.number,
                                   pairwise.interaction.number,
                                   saturated.number,
                                   baselineonly.second.cat,
                                   main.effect.second.cat,
                                   pairwise.interaction.second.cat,
                                   saturated.second.cat,
                                   M,
                                   total.covar.number,
                                   full.second.stage.names,
                                   covar.names,
                                   delta){

  max.number.second.stage.parameter <- 0
  if(baselineonly.number!=0){
    max.number.second.stage.parameter <- baselineonly.second.cat
  }
  if(main.effect.number!=0){
    max.number.second.stage.parameter <- main.effect.second.cat
  }
  if(pairwise.interaction.number!=0){
    max.number.second.stage.parameter <- pairwise.interaction.second.cat
  }
  if(saturated.number!=0){
    max.number.second.stage.parameter <- saturated.second.cat
  }


  second.stage.mat <- matrix(NA,max.number.second.stage.parameter,(total.covar.number-1))

  delta.no.inter <- delta[(M+1):length(delta)]
  ind.delta <- 0
  ind.covar <- 0
  if(baselineonly.number!=0){
    baselineonly.second.stage.delta <- matrix(NA,max.number.second.stage.parameter,baseline.number)
    baselineonly.second.stage.delta[1,(ind.covar+1):(ind.covar+baselineonly.number)] <- delta.no.inter[(ind.covar+1):(ind.covar+baselineonly.number)]
    second.stage.mat[,(ind.covar+1):(ind.covar+baselineonly.number)]<- baselineonly.second.stage.delta
    ind.covar <- ind.covar+baselineonly.number
    ind.delta <- ind.delta + ind.covar*baseline.number
  }

  if(main.effect.number!=0){

    for(i in 1:main.effect.number){
      print(i)
      ind.covar <- ind.covar+1
      second.stage.mat[1:main.effect.second.cat,ind.covar]<-
        delta.no.inter[(ind.delta+1):(ind.delta+main.effect.second.cat)]
      ind.delta <- ind.delta + main.effect.second.cat
    }
  }


  if(pairwise.interaction.number!=0){

    for(i in 1:pairwise.interaction.number){
      ind.covar <- ind.covar+1
      second.stage.mat[1:pairwise.interaction.second.cat,(ind.covar+1):(ind.covar+pairwise.interaction.number)]<-
        delta.no.inter[(ind.delta+1):(ind.delta+pairwise.interaction.second.cat)]
      ind.delta <- ind.delta + pairwise.interaction.second.cat
    }
  }

  if(saturated.number!=0){

    for(i in 1:saturated.number){
      ind.covar <- ind.covar+1
      second.stage.mat[1:saturated.second.cat,(ind.covar+1):(ind.covar+saturated.number)]<-
        delta.no.inter[(ind.delta+1):(ind.delta+saturated.second.cat)]
      ind.delta <- ind.delta + saturated.second.cat
    }
  }

  colnames(second.stage.mat) <- covar.names
  rownames(second.stage.mat) <- full.second.stage.names[1:max.number.second.stage.parameter]
  return(second.stage.mat)

}




