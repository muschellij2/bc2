#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
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
                          missingTumorIndicator = NULL){
  if(is.null(missingTumorIndicator)==1){
    return(Mvpoly(y,
                  baselineonly,
                  additive,
                  pairwise.interaction,
                  saturated))
  }else{
   return(EMmvpoly(y,
                               baselineonly,
                               additive,
                               pairwise.interaction,
                               saturated,
                               missingTumorIndicator))

  }
}
