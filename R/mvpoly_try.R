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



