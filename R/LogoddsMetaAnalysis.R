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
  info.meta <- info1+info2
  sigma.meta <- solve(info1+info2)
  logodds.meta <- sigma.meta%*%(info1%*%logodds1+info2%*%logodds2)

  return(list(logodds.meta = logodds.meta,
         info.meta = info.meta))

}
