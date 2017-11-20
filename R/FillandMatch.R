
#' Title
#'
#' @param reference
#' @param all
#'
#' @return
#' @export
#'
#' @examples
FillandMatch <- function(reference,all){
  idx.fil <- which(all%in%reference)
  all.fill <- all[idx.fil]
  idx.match <- match(reference,all.fill)
  return(list(idx.fill,idx.match))
}
