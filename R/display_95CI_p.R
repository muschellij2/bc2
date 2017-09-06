#' Title
#'
#' @param delta_odds
#' @param delta_test_low_odds
#' @param delta_test_high_odds
#' @param p_value
#'
#' @return
#' @export
#'
#' @examples
display_95CI_p = function(delta_odds,delta_test_low_odds,delta_test_high_odds,p_value){
  result = NULL
  for(i in 1:length(delta_odds)){
    result= c(result,paste0(delta_odds[i],"(",delta_test_low_odds[i],"-",
                            delta_test_high_odds[i],")"),
              p_value[i])
  }
  return(result)
}
