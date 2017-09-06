#' Title
#'
#' @param z
#'
#' @return
#' @export
#'
#' @examples
p_value_function <- function(z){
  result <- NULL
  for(i in 1:length(z)){

    result <- c(result,2*(pnorm(-abs(z[i]))))

  }
  return(result)
}
