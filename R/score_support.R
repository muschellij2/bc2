
#' Title
#'
#' @param pxx
#' @param x_covar
#' @param z_design
#' @param y_em
#'
#' @return
#' @export
#'
#' @examples
score_support <- function(pxx,x_covar,z_design,y_em){
  N <- nrow(x_covar)
  p_covar <- ncol(x_covar)

  ###need improvement
  M <- length(y_em)/N
  ret_p <- pxx
  # sof <- "getW.so"
  # dyn.load(sof)
  W <- rep(0,N*M*M)
  temp <- .C("Weighted_W", as.numeric(ret_p), as.numeric(W),as.integer(N), as.integer(M))
  W.c <- temp[[2]]
  temp_mis <- .C("Weighted_W", as.numeric(y_em), as.numeric(W),as.integer(N), as.integer(M))
  W.mis <- temp_mis[[2]]
  W <- W.c-W.mis
  x_covar <- cbind(1,x_covar)
  p_col <- ncol(x_covar)

  # sof1 <- "getXpWXp.so"
  # dyn.load(sof1)
  XpWXp <- rep(0,p_col^2*M^2)
  XpWXp_C <- .C("getXpWXp",as.numeric(W),as.numeric(x_covar),as.integer(N),as.integer(M),
                as.numeric(XpWXp),as.integer(p_col))
  XpWXp <- matrix(XpWXp_C[[5]],ncol(x_covar)*M,ncol(x_covar)*M)


  # sof <- "getWXp.so"
  # dyn.load(sof)

  ret <- rep(0,N*M*M*p_col)
  Wxp <- .C("getWXp",as.numeric(W),as.numeric(x_covar),as.integer(N),as.integer(M),as.numeric(ret),as.integer(p_col))
  WxpVec <- Wxp[[5]]



  z_covar <- ztozall.fun(p_covar,z_design,M)




  z_test <- z_design
  zpz_inverse <- solve(t(z_covar)%*%XpWXp%*%z_covar)

  return(list(z_covar=z_covar,z_test=z_test,zpz_inverse=zpz_inverse,W=W,WxpVec=WxpVec,ret_p=ret_p,
              x_covar = x_covar[,-1]))



  # score_try <- t(z_design)%*%XYP
  #
  #
  # infor_try <- t(z_design)%*%XWX%*%z_design-t(z_design)%*%XWXp_mat%*%z_covar%*%solve(t(z_covar)%*%XpWXp%*%z_covar)%*%t(z_covar)%*%t(XWXp_mat)%*%z_design
  #
  # score_test_value_try <- t(score_try)%*%solve(infor_try)%*%score_try

}
