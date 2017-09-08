#' Title
#'
#' @param y_em
#' @param p
#' @param missing.vec
#' @param missing.mat
#' @param nparm
#'
#' @return
#' @export
#'
#' @examples
LogLikelihoodwithAIC <- function(y_em,p,missing.vec,nparm){
  n <- nrow(y_em)
  M <- ncol(y_em)
  all.vec <- c(1:n)
  complete.vec <- which(all.vec%in%missing.vec==0)

  complete.loglikelihood <- ComputeLogLikelihood(y_em,p)
  AIC.complete.loglikelihood <- 2*nparm-2*complete.loglikelihood
  y.em.complete <- y_em[complete.vec,]

  p.mat <- matrix(p,nrow=n,ncol=M)
  p.mat.complete <- p.mat[complete.vec,]
  p.complete <- as.vector(p.mat.complete)
  loglikehood.for.complete.obs <- ComputeLogLikelihood(y.em.complete,p.complete)
  AIC.for.complete.obs <- 2*nparm - 2* loglikehood.for.complete.obs
  return(list(complete.loglikelihood=complete.loglikelihood,
              AIC.complete.loglikelihood = AIC.complete.loglikelihood,
              loglikehood.for.complete.obs=loglikehood.for.complete.obs,
              AIC.for.complete.obs=AIC.for.complete.obs))
}

#' Title
#'
#' @param y_em
#' @param p
#'
#' @return
#' @export
#'
#' @examples
ComputeLogLikelihood <- function(y_em,p){
  n <- nrow(y_em)
  M <- ncol(y_em)
  complete.loglikelihood <- 0
  for(i in 1:n){
    p.temp <- p[i+n*(0:(M-1))]
    y.temp <- y_em[i,]
    complete.loglikelihood <- complete.loglikelihood+
      crossprod(y.temp,log(p.temp))+(1-sum(y.temp))*log(1-sum(p.temp))
  }
  return(complete.loglikelihood)
}


