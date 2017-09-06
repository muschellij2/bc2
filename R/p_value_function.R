#' Title
#'
#' @param z
#'
#' @return
#' @export
#'
#' @examples
PvalueFunction <- function(z){
  result <- NULL
  for(i in 1:length(z)){

    result <- c(result,2*(pnorm(-abs(z[i]))))

  }
  return(result)
}
