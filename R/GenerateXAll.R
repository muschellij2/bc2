#' Title
#'
#' @param y
#' @param baselineonly
#' @param main.effect
#' @param pairwise.interaction
#' @param saturated
#'
#' @return
#' @export
#'
#' @examples
GenerateXAll <- function(y,baselineonly,main.effect,pairwise.interaction,saturated){
  n = nrow(y)
  ###initial x.all to use cbind
  x.all = rep(1,n)
  if(is.null(baselineonly)==0){
    x.all = cbind(x.all,baselineonly)
  }
  if(is.null(main.effect)==0){
    x.all = cbind(x.all,main.effect)
  }
  if(is.null(pairwise.interaction)==0){
    x.all = cbind(x.all,pairwise.interaction)
  }
  if(is.null(saturated)==0){
    x.all = cbind(x.all,saturated)
  }
  x.all <- x.all[,-1]
  return(x.all)
}

