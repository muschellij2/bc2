
#' Title
#'
#' @param y
#' @param x_test
#' @param score_support_result
#'
#' @return
#' @export
#'
#' @examples
score_test_mis <- function(y,x_test,score_support_result){
  if(nrow(y)<450){
    # if(is.null(x_covar)){
    #   yy <- as.vector(y)
    #   N <- nrow(x_test)
    #   p_covar <- 0
    #   M <- ncol(y)
    #   I_M <- diag(M)
    #   estimate <- Mvpoly(delta0,y,NULL,z_design)
    #   mu <- estimate$mu
    #   W <- estimate$W
    #
    #   # z_all <- NULL
    #   # for(i in 1:M){
    #   #   z_all <- rbind(z_all,kronecker(diag(p_covar+1),t(z[i,])))
    #   # }
    #   z_covar <- diag(M)
    #
    #   z_test <- z_design
    #   xx_test <- kronecker(I_M,x_test)
    #   xxz_test <- xx_test%*%z_test
    #   score <- t(xxz_test)%*%(yy-mu)
    #   xx_covar <- kronecker(I_M,rep(1,N))
    #   xxz_covar <- xx_covar%*%z_covar
    #   info <- t(xxz_test)%*%W%*%xxz_test-t(xxz_test)%*%W%*%xxz_covar%*%solve(infor_covar)%*%t(xxz_covar)%*%W%*%xxz_test
    #   score_test_result <- t(score)%*%solve(info)%*%score
    #   p_value <- 1-pchisq(as.numeric(score_test_result),df=ncol(z_test))
    #   return(p_value)
    #
    # }

    yy <- as.vector(y)
    N <- nrow(as.matrix(x_test))
    p_covar <- ncol(x_covar)
    M <- ncol(y)
    I_M <- diag(M)

    epi <- 1e-04
    #delta_old <- rep(0,length(delta0))
    delta_old <- delta0
    ##EM algorithm
    z_covar <- ztozall.fun(ncol(x_covar),z_design,M)



    for(niter in 1:100){

      y_em <- prob_fitting(delta_old,y_pheno,cbind(1,x_covar),z_standard,z_covar)
      estimate <-  Mvpoly(delta_old,y_em,x_covar,z_design)
      delta_new <-estimate[[1]]
      error_f(delta_old,delta_new)
      if(error_f(delta_old,delta_new)<epi){

        print("Converge!!!!!")
        break
      }
      delta_old <- delta_new
      print(niter)

    }




    #infor_mis_c <- infor_mis(y_em,x,z_design)
    mu <- estimate$mu
    W_c <- weighted_matrix(mu,N)
    W_mis <- weighted_matrix(y_em,N)
    W <- W_c-W_mis

    z_test <- z_design
    xx_test <- kronecker(I_M,x_test)
    xxz_test <- xx_test%*%z_test
    score <- t(xxz_test)%*%(yy-mu)
    x_covar <- cbind(1,x_covar)
    xx_covar <- kronecker(I_M,x_covar)
    xxz_covar <- xx_covar%*%z_covar
    infor_covar <- t(xxz_covar)%*%W%*%xxz_covar


    info <- t(xxz_test)%*%W%*%xxz_test-t(xxz_test)%*%W%*%xxz_covar%*%solve(infor_covar)%*%t(xxz_covar)%*%W%*%xxz_test
    score_test_result <- t(score)%*%solve(info)%*%score
    p_value <- 1-pchisq(as.numeric(score_test_result),df=ncol(z_test))
    return(p_value)
  }else{

    yy <- as.vector(y)
    z_covar <- score_support_result$z_covar
    z_test <- score_support_result$z_test
    x_covar <- score_support_result$x_covar
    M <- nrow(z_test)
    N <- length(y)/M
    zpz_inverse <- score_support_result$zpz_inverse
    W <- score_support_result$W
    WxpVec <- score_support_result$WxpVec
    ret_p <- score_support_result$ret_p
    z_design <- z_test

    x_covar <- cbind(1,x_covar)
    p_col <- ncol(x_covar)

    sof <- "getXWXp.so"
    dyn.load(sof)
    ret <- rep(0,M^2*p_col)
    XWXp <- .C("getXWXp",as.numeric(WxpVec),as.numeric(x_test),as.integer(N),as.integer(M),as.numeric(ret),as.integer(p_col))
    XWXpVec <- XWXp[[5]]
    XWXp_mat_temp <- matrix(XWXpVec,M)
    XWXp_mat <- XWXp_mat_temp
    for(ind in 1:p_col){
      XWXp_mat[,seq(ind,M*p_col,p_col)] <- XWXp_mat_temp[,(ind-1)*M+(1:M)]
    }
    sof <- "getXWX.so"
    dyn.load(sof)
    ret = rnorm(M^2)
    XWX_c <- .C("getXWX",as.numeric(W),as.numeric(x_test),as.integer(N),as.integer(M),
                as.numeric(ret))
    XWX <- matrix(XWX_c[[5]],M,M)
    sof2 <- "getXY.so"
    dyn.load(sof2)
    XY <- rep(0,M)
    pxx <- yy-ret_p
    XY_c <- .C("getXY",as.numeric(x_test),as.numeric(pxx),as.numeric(XY),as.integer(N),as.integer(M))
    XYP <- XY_c[[3]]


    score_c <- t(z_design)%*%XYP


    infor_c <- t(z_design)%*%XWX%*%z_design-t(z_design)%*%XWXp_mat%*%z_covar%*%zpz_inverse%*%t(z_covar)%*%t(XWXp_mat)%*%z_design

    score_test_value_c <- t(score_c)%*%solve(infor_c)%*%score_c

    p_value <- 1-pchisq(as.numeric(score_test_value_c),df=ncol(z_test))

    return(list(score_c=score_c,infor_c=infor_c,p_value=p_value))


  }
}
