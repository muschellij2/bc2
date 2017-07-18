#This is a function implement in Rcpp to make it faster to calculate the weighted matrix during the iteration

#This Mvpoly needs four parameter to start with:

#First is delta0,

#Second is the  For example, if there are three characters for a specific cancer; HER+-;PR+-;ER+-; So there are 8 possible categories; In this case, M=8; If the first individual doesn't have disease, then he is (0,0,0,0,0,0,0,0); If the second individual is case with HER+,PR-,ER-, then he is (1,0,0,0,0,0,0,0)

#The third is x

#The fourth is z_design matrix.
#' Title
#'
#' @param delta0 this is the starting value for all the second stage parameters.
#' You need to assign a starting value to make the iteration algorithm begin
#' @param y disease status matrix y; This is an N by M matrix; N is the total number of samples; M is the total number of possible categories;
#' @param x matrix, which contains all the covariates you want to fit into the model;
#' @param z_design This matrix maps the relationship between the first stage parameters beta
#' and second stage parameters theta for a specific covariate.
#'  Notice:For a specific covariate! You don't need to write the matrix between
#'  first stage paramters beta and second stage paramters theta for all the covariate.
#'  The function will generate everything for you. The function use the same design matrix for every covariate in the model;For example, for the HER+-,ER+-,PR+- case, the z matrix is: the z_design matrix in simulation1.R
#'  #This is a shortcoming of this version of code. In the furture,
#'  I could develop the version that user could specific different z matrix for different covariates.
#'
#' @return
#' @export
#'
#' @examples
Mvpoly <- function(delta0, y, x, z_design) {
  if (nrow(y) < 450) {
    if (is.null(x)) {
      N <- nrow(y)
      x <- as.matrix(rep(1, N))
      p <- ncol(x)
      M <- ncol(y)
      I_M <- diag(M)
      I_N <- diag(N)
      yy <- as.vector(y)
      z <- diag(M)
      #xx <- kronecker(I_M,x)
      xx <- kronecker(I_M, x)
      one_M <- rep(1, M)
      #oneM_x_IN <- kronecker(one_M,I_N)
      xxz <- xx %*% z
      niter <- 1
      rerror <- 1
      delta <- delta0
      tol <- 1e-06
      #A <- matrix(0,nrow=N*M,ncol=N)
      #A2 <- matrix(0,nrow=N*M,ncol=N*M)
      #A <- Matrix(0,nrow = N*M,ncol=N,sparse=T)
      #A2 <- Matrix(0,nrow=N*M,ncol=N*M,sparse = T)
      yy <- as.vector(y)
      while ((niter < 100) & (rerror > tol)) {
        beta <- z %*% delta
        lxx <- xx %*% beta
        pxx <- exp(lxx)
        pxx <- matrix(pxx, N, M)
        spxx <- rowSums(pxx)
        pxx <- sweep(pxx, 1, (1 + spxx), "/")
        pxx <- as.vector(pxx)
        #A <- D%*% oneM_x_IN
        # for(i in 1:M){
        #   A[(1+(i-1)*N):(N*i),] <- D[(1+(i-1)*N):(N*i),(1+(i-1)*N):(N*i)]
        # }

        #A2 <- crossprod(t(A),t(A))
        # for(i in 1:M){
        #   for(j in 1:M){
        #     A2[(1+(i-1)*N):(N*i),(1+(j-1)*N):(N*j)] <-
        #       diag_times(A[(1+(i-1)*N):(N*i),],
        #                  A[(1+(j-1)*N):(N*j),])
        #
        #   }
        # }
        W <- weighted_matrix(pxx, N)

        W_y <- yy - pxx + W %*% lxx
        Infor_M <- t(xxz) %*% W %*% xxz
        delta <- solve(Infor_M, t(xxz) %*% W_y)
        rerror <- max(abs(delta - delta0) / max(abs(delta0), 0.1))
        niter <- niter + 1
        delta0 <- as.vector(delta)
        print(delta0)
        print(niter)
        if (any(is.nan(delta0)) == T) {
          print("the algorithm doesn't converge")
          break
        }
      }
      # score <- score_f(xxz,y,pxx)
      result <- list(delta = matrix(delta, ncol = p),
                     mu = pxx,
                     W = W)
      return(result)
    }
    N <- nrow(x)
    p_col <- ncol(x)
    x <- cbind(1, x)
    p <- ncol(x)
    M <- ncol(y)
    I_M <- diag(M)
    I_N <- diag(N)
    yy <- as.vector(y)
    z_all_temp <- NULL

    for (i in 1:M) {
      z_all_temp <-
        rbind(z_all_temp, kronecker(diag(p_col), t(z_design[i, ])))
    }

    z_all <-
      matrix(0,
             nrow = M * (p_col + 1),
             ncol = M + p_col * ncol(z_design))
    for (i in 1:M) {
      z_all[1 + (i - 1) * (p_col + 1), i] <- 1
    }
    for (i in 1:M) {
      z_all[(2 + (i - 1) * (p_col + 1)):(i * (p_col + 1)), (M + 1):ncol(z_all)] <-
        z_all_temp[(1 + (i - 1) * p_col):(i * p_col), ]
    }


    z <- z_all


    #xx <- kronecker(I_M,x)
    xx <- kronecker(I_M, x)
    one_M <- rep(1, M)
    #oneM_x_IN <- kronecker(one_M,I_N)
    xxz <- xx %*% z
    niter <- 1
    rerror <- 1
    delta <- delta0
    tol <- 1e-06
    #A <- matrix(0,nrow=N*M,ncol=N)
    #A2 <- matrix(0,nrow=N*M,ncol=N*M)
    #A <- Matrix(0,nrow = N*M,ncol=N,sparse=T)
    #A2 <- Matrix(0,nrow=N*M,ncol=N*M,sparse = T)
    #yy <- as.vector(y)
    while ((niter < 100) & (rerror > tol)) {
      beta <- z %*% delta
      lxx <- xx %*% beta
      pxx <- exp(lxx)
      pxx <- matrix(pxx, N, M)
      spxx <- rowSums(pxx)
      pxx <- sweep(pxx, 1, (1 + spxx), "/")
      pxx <- as.vector(pxx)
      #A <- D%*% oneM_x_IN
      # for(i in 1:M){
      #   A[(1+(i-1)*N):(N*i),] <- D[(1+(i-1)*N):(N*i),(1+(i-1)*N):(N*i)]
      # }

      #A2 <- crossprod(t(A),t(A))
      # for(i in 1:M){
      #   for(j in 1:M){
      #     A2[(1+(i-1)*N):(N*i),(1+(j-1)*N):(N*j)] <-
      #       diag_times(A[(1+(i-1)*N):(N*i),],
      #                  A[(1+(j-1)*N):(N*j),])
      #
      #   }
      # }
      W <- weighted_matrix(pxx, N)

      W_y <- yy - pxx + W %*% lxx
      Infor_M <- t(xxz) %*% W %*% xxz
      print(Infor_M)
      delta <- solve(Infor_M, t(xxz) %*% W_y)
      rerror <- max(abs(delta - delta0) / max(abs(delta0), 0.1))
      niter <- niter + 1
      delta0 <- as.vector(delta)
      print(delta0)
      print(niter)
      if (any(is.nan(delta0)) == T) {
        print("the algorithm doesn't converge")
        break
      }
    }
    # score <- score_f(xxz,y,pxx)
    result <- list(
      delta = as.matrix(delta),
      mu = pxx,
      W = W,
      infor = Infor_M
    )
    return(result)
  }


  # sof <- "source1.so"
  # sof = load_so(sof)
  # dyn.load(sof)

  N <- nrow(x)
  #x <- cbind(1,x)
  #p <- ncol(x)
  M <- ncol(y)
  NCAT   <- ncol(z_design) - 1
  NCATP1 <- NCAT + 1
  NCOV   <- ncol(x)
  NM     <- N * M
  tol    <- 1e-6
  nparm  <-  NCATP1 * (NCOV) + M
  deltai <- delta0

  NITER  <- 500
  Y <- as.vector(y)
  X <- as.vector(x)
  Z <- as.vector(as.matrix(z_design))



  debug     <- 1
  ret_rc    <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999, nparm ^ 2))
  ret_p <- as.numeric(rep(0, NM))

  temp <-
    .C(
      "Mvpoly",
      as.numeric(deltai),
      as.integer(nparm),
      as.numeric(Y),
      as.numeric(X),
      as.numeric(Z),
      as.integer(N),
      as.integer(M),
      as.integer(NCAT),
      as.integer(NCOV),
      as.integer(NITER),
      as.numeric(tol),
      as.integer(debug),
      ret_rc = ret_rc,
      ret_delta = ret_delta,
      ret_info = ret_info,
      ret_p = ret_p,
      PACKAGE = "bc2"
    )
  info <- matrix(unlist(temp$ret_info), nparm, nparm)
  result <- list(temp$ret_delta, info,
                 temp$ret_p)
  return(result)

}















# sof <- "source.so"
# dyn.load(sof)
#
#
#
# nparm <- length(delta0)
# NITER <- 100
# tol   <- 1e-4
# NCOV  <- ncol(x_covar)
# ncat  <- NCOL(z_design)-1
# DEBUG     <- 2
# ret_rc    <- as.integer(1)
# ret_delta <- as.numeric(rep(-9999, nparm))
#
# ret_p     <- as.numeric(1)
# M <- ncol(y)
# n <- nrow(y)
#
# temp <- .C("score_test", as.numeric(delta0), as.integer(nparm), as.numeric(as.vector(y)),
#            as.numeric(x_covar), as.numeric(x_test), as.numeric(z_design), as.integer(n), as.integer(M),
#            as.integer(ncat), as.integer(NCOV), as.integer(NITER), as.numeric(tol),
#            as.integer(DEBUG), ret_rc=ret_rc, ret_delta=ret_delta, ret_p=ret_p)
# return(temp$ret_p)
# }

#}


# diag_times <- function(A,B){
#   n <- nrow(A)
#   result <- matrix(0,nrow=n,ncol=n)
#   for(i in 1:n){
#     result[i,i] <- A[i,i]*B[i,i]
#   }
#   return(result)
# }

# diag_times <- function(A,B){
#   n <- nrow(A)
#   result <- matrix(0,nrow=n,ncol=n)
#   for(i in 1:n){
#     result[i,i] <- A[i,i]*B[i,i]
#   }
#   return(result)
# }
