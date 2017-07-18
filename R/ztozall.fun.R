


#' Title
#'
#' @param p_covar
#' @param z_design
#' @param M
#'
#' @return
#' @export
#'
#' @examples
ztozall.fun <- function(p_covar,z_design,M){

  z_covar_temp <- NULL

  for(i in 1:M){
    z_covar_temp <- rbind(z_covar_temp,kronecker(diag(p_covar),t(z_design[i,])))
  }

  z_covar <- matrix(0,nrow = M*(p_covar+1),ncol= M+p_covar*ncol(z_design))
  for(i in 1:M){
    z_covar[1+(i-1)*(p_covar+1),i] <- 1
  }
  for(i in 1:M){
    z_covar[(2+(i-1)*(p_covar+1)):(i*(p_covar+1)),(M+1):ncol(z_covar)] <-
      z_covar_temp[(1+(i-1)*p_covar):(i*p_covar),]
  }
  return(z_covar)

}
