
infor_mis <- function(y,x,z){


  sof <- "source_mis.so"
  dyn.load(sof)

  N <- nrow(x)
  #x <- cbind(1,x)
  #p <- ncol(x)
  M <- ncol(y)
  NCAT   <- ncol(z_design)-1
  NCATP1 <- NCAT + 1
  NCOV   <- ncol(x)
  NM     <- N*M
  tol    <- 1e-6
  nparm  <-  NCATP1*(NCOV)+M

  NITER  <- 500
  Y <- as.vector(y)
  X <- as.vector(x)
  Z <- as.vector(as.matrix(z_design))



  debug     <- 1
  ret_rc    <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999,nparm^2))

  temp <- .C("mis_infor", as.integer(nparm), as.numeric(Y), as.numeric(X), as.numeric(Z),
             as.integer(N), as.integer(M), as.integer(NCAT), as.integer(NCOV), as.integer(NITER),
             as.numeric(tol), as.integer(debug), ret_rc=ret_rc,
             ret_info = ret_info)
  info <- matrix(unlist(temp$ret_info),nparm,nparm)
  return(info)

}
