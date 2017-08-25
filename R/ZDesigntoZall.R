#' Title
#'
#' @param baselineonly
#' @param main.effect
#' @param pairwise.interaction
#' @param saturated
#' @param z.design.baselineonly
#' @param z.design.main.effect
#' @param z.design.pairwise.interaction
#' @param z.design.saturated
#'
#' @return
#' @export
#'
#' @examples
ZDesigntoZall <- function(baselineonly,
                          main.effect,
                          pairwise.interaction,
                          saturated,
                          z.design.baselineonly,
                          z.design.main.effect,
                          z.design.pairwise.interaction,
                          z.design.saturated) {
  M <- nrow(z.design.main.effect)
  ##the number of covariates in different potential model structures
  baselineonly.number <- CountCovarNumber(baselineonly)
  main.effect.number <- CountCovarNumber(main.effect)
  pairwise.interaction.number <- CountCovarNumber(pairwise.interaction)
  saturated.number <- CountCovarNumber(saturated)
  ###second.stage.category for different model structures
  baselineonly.second.cat <- 1
  main.effect.second.cat <- ncol(z.design.main.effect)
  pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
  saturated.second.cat <- ncol(z.design.saturated)
  ###1 for intercept
  total.covar.number <- 1+ baselineonly.number+main.effect.number+
    pairwise.interaction.number+saturated.number

  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+baselineonly.number*baselineonly.second.cat+
                                                          main.effect.second.cat*main.effect.number+
                                                          pairwise.interaction.second.cat*pairwise.interaction.number)+saturated.second.cat*saturated.number)

  for(i in c("intercept","baselineonly",
             "main.effect","pairwise.interection",
             "satuared")){
    ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
    if(i=="intercept"){
      ###row start and column start point for this category
      row.start <- 0
      column.start <- 0
      for(j in 1:M){
        z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
      }
    }else if(i=="baselineonly"){
      column.start = M
      row.start <- 1
      ###test whether there is any baselineonly variable
      if(baselineonly.number!=0){
        for(j in 1:M){
          for(k in 1:baselineonly.number){
            z.all[row.start+k+(j-1)*total.covar.number,(column.start+k)] <- 1
          }
        }
      }
    }else if(i=="main.effect"){
      column.start <- M+baselineonly.number
      row.start <- 1+baselineonly.number
      if(main.effect.number!=0){
        for(j in 1:M){
          for(k in 1:main.effect.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*main.effect.second.cat+1):
                    (column.start+k*main.effect.second.cat)] <- as.matrix(z.design.main.effect[j,])
          }
        }
      }
    }else if(i == "pairwise.interaction"){
      column.start <- M+baselineonly.number+main.effect.number*main.effect.second.cat
      row.start <- 1+baselineonly.number+main.effect.number
      if(pairwise.interaction.number!=0){
        for(j in 1:M){
          for(k in 1:pairwise.interaction.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*pairwise.interaction.second.cat+1):
                    (column.start+k*pairwise.interaction.second.cat)] <- as.matrix(z.design.pairwise.interaction[j,])
          }
        }
      }
    }else {
      column.start <- M+baselineonly.number+main.effect.number*main.effect.second.cat+
        pairwise.interaction.number*pairwise.interaction.second.cat
      row.start <- 1+baselineonly.number+main.effect.number+pairwise.interaction.number
      if(saturated.number!=0){
        for(j in 1:M){
          for(k in 1:saturated.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*saturated.second.cat+1):
                    (column.start+k*saturated.second.cat)] <- as.matrix(z.design.saturated[j,])
          }
        }
      }
    }



  }
  return(z.all)
}
