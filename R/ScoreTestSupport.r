
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
ScoreTestSupport <- function(y,
                               baselineonly=NULL,
                               additive=NULL,
                               pairwise.interaction=NULL,
                               saturated=NULL,
                               missingTumorIndicator = NULL,
                             missingDataAlgorithm = "EM"
                             ){
 if(is.null(missingTumorIndicator)==1){
  return(CompleteCasesScoreTestSupport(y,
                                                   baselineonly,
                                                   additive,
                                                   pairwise.interaction,
                                                   saturated))
 }else{
   if(missingDataAlgorithm=="EM"){
    return(EMScoreTestSupport(y,
                                          baselineonly,
                                          additive,
                                          pairwise.interaction,
                                          saturated,
                                          missingTumorIndicator))
   }else if(missingDataAlgorithm=="OneStepMLE"){

   }
 }

}
