#' Title
#'
#' @param z.design
#' @param z.design.additive
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ZSelfDesigntoZall <- function(z.design,
                              z.design.additive,
                              x){

  M <- nrow(z.design.additive)
  if(is.vector(x)==1){
    x = matrix(x,ncol = 1)
  }
  total.covar.number <- ncol(x)+1
  total.covar.number.no.inter <- ncol(x)
  second.stage.cat <- ncol(z.design)
  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+total.covar.number.no.inter*second.stage.cat))

  ##setup intercept
  ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept

  row.start <- 0
  column.start <- 0
  for(j in 1:M){
    z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
  }
  ##set up all other covariates
      column.start = M
      row.start <- 1
      ###test whether there is any baselineonly variable

      for(j in 1:M){
        for(k in 1:(total.covar.number-1)){
          z.all[row.start+k+(j-1)*total.covar.number,
                (column.start+(k-1)*second.stage.cat+1):
                  (column.start+k*second.stage.cat)] <- as.matrix(z.design[j,])
        }
      }
       return(z.all)
}
