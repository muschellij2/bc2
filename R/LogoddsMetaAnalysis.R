#' Title
#'
#' @param logodds1
#' @param info1
#' @param logodds2
#' @param info2
#'
#' @return
#' @export
#'
#' @examples
LogoddsMetaAnalysis <- function(logodds1,info1,logodds2,info2){
  w1 <- solve(info1)
  w2 <- solve(info2)
  logodds.meta <- solve(w1+w2)%*%(w1%*%logodds1+w2%*%logodds2)
  info.meta <- solve(w1+w2)

  return(list(logodds.meta = logodds.meta,
         info.meta = info.meta))

}
