#' Title
#'
#' @param all
#' @param fine_mapping
#'
#' @return
#' @export
#'
#' @examples
get_fine_mapping_id<- function(all,fine_mapping){
CHR.all <- all$CHR
position.all <- all$position
for(i in 1:nrow(fine_mapping)){
  #print(i)
  chr_temp <- CHR[i]
  start_temp <- start[i]
  end_temp <- end[i]
  idx <- which(CHR.all==chr_temp&position.all>=start_temp&
                 position.all<=end_temp)
  temp.known.flag <- rep(i,length(idx))
  idx_cut <- c(idx_cut,idx)
  known.flag <- c(known.flag,temp.known.flag)
}
return(list(idx_cut,known.flag))

}
