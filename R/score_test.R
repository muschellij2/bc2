

###This score_test contains 5 paramters;
#delta0,y and z is the same as the mvpoly funciton; here z means the z_design matrix;
#x_test is the covariate vector you want to test in your model. Your null hypotheis is all the second stage parameters theta are 0 for this covariate
#x_covar is the matrix for all the other covariates you want to fit into the model;
#' Title
#'
#' @param delta0
#' @param y
#' @param x_test
#' @param x_covar
#' @param z_design
#'
#' @return
#' @export
#'
#' @examples
score_test <- function(delta0, y, x_test, x_covar, z_design) {
  if (nrow(y) < 495) {
    if (is.null(x_covar)) {
      yy <- as.vector(y)
      N <- nrow(x_test)
      p_covar <- 0
      M <- ncol(y)
      I_M <- diag(M)
      estimate <- Mvpoly(delta0, y, NULL, z_design)
      mu <- estimate$mu
      W <- estimate$W

      # z_all <- NULL
      # for(i in 1:M){
      #   z_all <- rbind(z_all,kronecker(diag(p_covar+1),t(z[i,])))
      # }
      z_covar <- diag(M)

      z_test <- z_design
      xx_test <- kronecker(I_M, x_test)
      xxz_test <- xx_test %*% z_test
      score <- t(xxz_test) %*% (yy - mu)
      xx_covar <- kronecker(I_M, rep(1, N))
      xxz_covar <- xx_covar %*% z_covar
      info <-
        t(xxz_test) %*% W %*% xxz_test - t(xxz_test) %*% W %*% xxz_covar %*% solve(infor_covar) %*%
        t(xxz_covar) %*% W %*% xxz_test
      score_test_result <- t(score) %*% solve(info) %*% score
      p_value <-
        1 - pchisq(as.numeric(score_test_result), df = ncol(z_test))
      return(p_value)

    }
    yy <- as.vector(y)
    N <- nrow(x_test)
    p_covar <- ncol(x_covar)
    M <- ncol(y)
    I_M <- diag(M)
    estimate <- Mvpoly(delta0, y, x_covar, z_design)
    mu <- estimate$mu
    W <- estimate$W
    infor_covar <- estimate$infor

    z_covar_temp <- NULL

    for (i in 1:M) {
      z_covar_temp <-
        rbind(z_covar_temp, kronecker(diag(p_covar), t(z_design[i, ])))
    }

    z_covar <-
      matrix(0,
             nrow = M * (p_covar + 1),
             ncol = M + p_covar * ncol(z_design))
    for (i in 1:M) {
      z_covar[1 + (i - 1) * (p_covar + 1), i] <- 1
    }
    for (i in 1:M) {
      z_covar[(2 + (i - 1) * (p_covar + 1)):(i * (p_covar + 1)), (M + 1):ncol(z_covar)] <-
        z_covar_temp[(1 + (i - 1) * p_covar):(i * p_covar), ]
    }




    z_test <- z_design
    xx_test <- kronecker(I_M, x_test)
    xxz_test <- xx_test %*% z_test
    score <- t(xxz_test) %*% (yy - mu)
    x_covar <- cbind(1, x_covar)
    xx_covar <- kronecker(I_M, x_covar)
    xxz_covar <- xx_covar %*% z_covar

    info <-
      t(xxz_test) %*% W %*% xxz_test - t(xxz_test) %*% W %*% xxz_covar %*% solve(infor_covar) %*%
      t(xxz_covar) %*% W %*% xxz_test
    score_test_result <- t(score) %*% solve(info) %*% score
    p_value <-
      1 - pchisq(as.numeric(score_test_result), df = ncol(z_test))
    return(p_value)
  } else{
    sof <- "source.so"
    dyn.load(sof)



    nparm <- length(delta0)
    NITER <- 100
    tol   <- 1e-4
    NCOV  <- ncol(x_covar)
    ncat  <- NCOL(z_design) - 1
    DEBUG     <- 2
    ret_rc    <- as.integer(1)
    ret_delta <- as.numeric(rep(-9999, nparm))

    ret_p     <- as.numeric(1)
    M <- ncol(y)
    n <- nrow(y)

    temp <-
      .C(
        "score_test",
        as.numeric(delta0),
        as.integer(nparm),
        as.numeric(as.vector(y)),
        as.numeric(x_covar),
        as.numeric(x_test),
        as.numeric(z_design),
        as.integer(n),
        as.integer(M),
        as.integer(ncat),
        as.integer(NCOV),
        as.integer(NITER),
        as.numeric(tol),
        as.integer(DEBUG),
        ret_rc = ret_rc,
        ret_delta = ret_delta,
        ret_p = ret_p
      )
    return(temp$ret_p)
  }

}
