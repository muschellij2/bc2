
#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
#' @param missingDataAlgorithm
#'
#' @return
#' @export
#'
#' @examples
TwoStageModel <- function(y,
                          baselineonly=NULL,
                          additive=NULL,
                          pairwise.interaction=NULL,
                          saturated=NULL,
                          missingTumorIndicator = NULL,
                          missingDataAlgorithm = "EM"){
  if(is.null(missingTumorIndicator)==1){
    return(Mvpoly(y,
                  baselineonly,
                  additive,
                  pairwise.interaction,
                  saturated))
  }else{
    if(missingDataAlgorithm=="EM"){
      return(EMmvpoly(y,
                      baselineonly,
                      additive,
                      pairwise.interaction,
                      saturated,
                      missingTumorIndicator))
    }else if(missingDataAlgorithm=="OneStepMLE"){
      return(OneStepMLE(y,
                                    baselineonly,
                                    additive,
                                    pairwise.interaction,
                                    saturated,
                                    missingTumorIndicator))
    }else{
      print(paste0("no missing data algorithm called ",missingDataAlgorithm))
    }


  }
}
