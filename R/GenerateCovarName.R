
#' Title
#'
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#'
#' @return
#' @export
#'
#' @examples
GenerateCovarName <- function(baselineonly,
                              additive,
                              pairwise.interaction,
                              saturated){
  full.names <- NULL
  if(is.null(baselineonly)==0){
    full.names<- c(full.names,colnames(baselineonly))
  }
  if(is.null(additive)==0){
    full.names<- c(full.names,colnames(additive))
  }
  if(is.null(pairwise.interaction)==0){
    full.names<- c(full.names,colnames(pairwise.interaction))
  }
  if(is.null(saturated)==0){
    full.names<- c(full.names,colnames(saturated))
  }
  return(full.names)
}
