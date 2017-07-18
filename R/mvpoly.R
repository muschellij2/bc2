library(Matrix)
library(Rcpp)
logit_inver <- function(x){
  return(exp(x)/(1+exp(x)))
}



####JUL 26 two stage modeling on 2004 JASA ########
#This code include two main function:Mvpoly and score_test function;
#Mvpoly gives the MLE under the two stage model
#Score_test gives the score test result under the NULL hypothesis that every second parameter theta for that specific test covariate is 0;


#This is a function implement in Rcpp to make it faster to calculate the weighted matrix during the iteration
cppFunction('NumericMatrix weighted_matrix(NumericVector p,int N){
            int NM=p.size();
            int M = NM/N;
            NumericMatrix w(NM,NM);
            for(int i=0; i < M; ++i){
            for(int j=0; j< M; ++j){
            if(i==j){
            for(int k=0; k<N;++k){
            w(N*i+k,N*i+k)=p(N*i+k)-p(N*i+k)*p(N*i+k);
            }
            }
            else{
            for(int k=0; k<N;++k){
            w(i*N+k,j*N+k)=
            -p(i*N+k)*p(j*N+k);
            }
            }
            }
            }
            return(w);
            }')

#This Mvpoly needs four parameter to start with:

#First is delta0, this is the starting value for all the second stage parameters. You need to assign a starting value to make the iteration algorithm begin

#Second is the disease status matrix y; This is an N by M matrix; N is the total number of samples; M is the total number of possible categories; For example, if there are three characters for a specific cancer; HER+-;PR+-;ER+-; So there are 8 possible categories; In this case, M=8; If the first individual doesn't have disease, then he is (0,0,0,0,0,0,0,0); If the second individual is case with HER+,PR-,ER-, then he is (1,0,0,0,0,0,0,0)

#The third is x matrix, which contains all the covariates you want to fit into the model;

#The fourth is z_design matrix. This matrix maps the relationship between the first stage parameters beta and second stage parameters theta for a specific covariate. Notice:For a specific covariate! You don't need to write the matrix between first stage paramters beta and second stage paramters theta for all the covariate. The function will generate everything for you. The function use the same design matrix for every covariate in the model;For example, for the HER+-,ER+-,PR+- case, the z matrix is: the z_design matrix in simulation1.R
#This is a shortcoming of this version of code. In the furture, I could develop the version that user could specific different z matrix for different covariates.
Mvpoly <- function(delta0,y,x,z_design){
  if(nrow(y)<450){
    if(is.null(x)){
      N <- nrow(y)
      x <- as.matrix(rep(1,N))
      p <- ncol(x)
      M <- ncol(y)
      I_M <- diag(M)
      I_N <- diag(N)
      yy <- as.vector(y)
      z <- diag(M)
      #xx <- kronecker(I_M,x)
      xx <- kronecker(I_M,x)
      one_M <- rep(1,M)
      #oneM_x_IN <- kronecker(one_M,I_N)
      xxz <- xx%*%z
      niter <- 1
      rerror <- 1
      delta <- delta0
      tol <- 1e-06
      #A <- matrix(0,nrow=N*M,ncol=N)
      #A2 <- matrix(0,nrow=N*M,ncol=N*M)
      #A <- Matrix(0,nrow = N*M,ncol=N,sparse=T)
      #A2 <- Matrix(0,nrow=N*M,ncol=N*M,sparse = T)
      yy <- as.vector(y)
      while((niter<100)&(rerror>tol)){
        beta <- z%*%delta
        lxx <- xx%*%beta
        pxx <- exp(lxx)
        pxx <- matrix(pxx,N,M)
        spxx <- rowSums(pxx)
        pxx <- sweep(pxx,1,(1+spxx),"/")
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
        W <- weighted_matrix(pxx,N)
        
        W_y <- yy-pxx+W%*%lxx
        Infor_M <- t(xxz)%*%W%*%xxz
        delta <- solve(Infor_M,t(xxz)%*%W_y)
        rerror <- max(abs(delta-delta0)/max(abs(delta0),0.1))
        niter <- niter+1
        delta0 <- as.vector(delta)
        print(delta0)
        print(niter)
        if(any(is.nan(delta0))==T){
          print("the algorithm doesn't converge")
          break
        }
      }
      # score <- score_f(xxz,y,pxx)
      result <- list(delta=matrix(delta,ncol=p),mu=pxx,W=W)
      return(result)
    }
    N <- nrow(x)
    p_col <- ncol(x)
    x <- cbind(1,x)
    p <- ncol(x)
    M <- ncol(y)
    I_M <- diag(M)
    I_N <- diag(N)
    yy <- as.vector(y)
    z_all_temp <- NULL
    
    for(i in 1:M){
      z_all_temp <- rbind(z_all_temp,kronecker(diag(p_col),t(z_design[i,])))
    }
    
    z_all <- matrix(0,nrow = M*(p_col+1),ncol= M+p_col*ncol(z_design))
    for(i in 1:M){
      z_all[1+(i-1)*(p_col+1),i] <- 1
    }
    for(i in 1:M){
      z_all[(2+(i-1)*(p_col+1)):(i*(p_col+1)),(M+1):ncol(z_all)] <- 
        z_all_temp[(1+(i-1)*p_col):(i*p_col),]
    }
    
    
    z <- z_all
    
    
    #xx <- kronecker(I_M,x)
    xx <- kronecker(I_M,x)
    one_M <- rep(1,M)
    #oneM_x_IN <- kronecker(one_M,I_N)
    xxz <- xx%*%z
    niter <- 1
    rerror <- 1
    delta <- delta0
    tol <- 1e-06
    #A <- matrix(0,nrow=N*M,ncol=N)
    #A2 <- matrix(0,nrow=N*M,ncol=N*M)
    #A <- Matrix(0,nrow = N*M,ncol=N,sparse=T)
    #A2 <- Matrix(0,nrow=N*M,ncol=N*M,sparse = T)
    #yy <- as.vector(y)
    while((niter<100)&(rerror>tol)){
      beta <- z%*%delta
      lxx <- xx%*%beta
      pxx <- exp(lxx)
      pxx <- matrix(pxx,N,M)
      spxx <- rowSums(pxx)
      pxx <- sweep(pxx,1,(1+spxx),"/")
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
      W <- weighted_matrix(pxx,N)
      
      W_y <- yy-pxx+W%*%lxx
      Infor_M <- t(xxz)%*%W%*%xxz
      print(Infor_M)
      delta <- solve(Infor_M,t(xxz)%*%W_y)
      rerror <- max(abs(delta-delta0)/max(abs(delta0),0.1))
      niter <- niter+1
      delta0 <- as.vector(delta)
      print(delta0)
      print(niter)
      if(any(is.nan(delta0))==T){
        print("the algorithm doesn't converge")
        break
      }
    }
    # score <- score_f(xxz,y,pxx)
    result <- list(delta=as.matrix(delta),mu=pxx,W=W,infor=Infor_M)
    return(result)
  }
  
  
  sof <- "source1.so"
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
  deltai <- delta0
  
  NITER  <- 500
  Y <- as.vector(y)
  X <- as.vector(x)
  Z <- as.vector(as.matrix(z_design))
  
  
  
  debug     <- 1
  ret_rc    <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999,nparm^2))
  ret_p <- as.numeric(rep(0,NM))
  
  temp <- .C("Mvpoly", as.numeric(deltai), as.integer(nparm), as.numeric(Y), as.numeric(X), as.numeric(Z), 
             as.integer(N), as.integer(M), as.integer(NCAT), as.integer(NCOV), as.integer(NITER), 
             as.numeric(tol), as.integer(debug), ret_rc=ret_rc, ret_delta=ret_delta,
             ret_info = ret_info,ret_p=ret_p)
  info <- matrix(unlist(temp$ret_info),nparm,nparm)
  result <- list(temp$ret_delta,info,
                 temp$ret_p)
  return(result)
  
}

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

getW <- function(p.complete,p.mis,N,M){
  sof <- "getW.so"
  dyn.load(sof)
  W.complete <- rep(0,N*M*M)
  temp <- .C("Weighted_W", as.numeric(p.complete), as.numeric(W.complete),as.integer(N), as.integer(M))
  W.complete <- temp[[2]]
  temp <- .C("Weighted_W", as.numeric(p.mis), as.numeric(W.complete),as.integer(N), as.integer(M))
  W.mis <- temp[[2]]
  return(W.complete-W.mis)
}


###This score_test contains 5 paramters;
#delta0,y and z is the same as the mvpoly funciton; here z means the z_design matrix;
#x_test is the covariate vector you want to test in your model. Your null hypotheis is all the second stage parameters theta are 0 for this covariate
#x_covar is the matrix for all the other covariates you want to fit into the model;
score_test <- function(delta0,y,x_test,x_covar,z_design){
  if(nrow(y)<495){
    if(is.null(x_covar)){
      yy <- as.vector(y)
      N <- nrow(x_test)
      p_covar <- 0
      M <- ncol(y)
      I_M <- diag(M)
      estimate <- Mvpoly(delta0,y,NULL,z_design)
      mu <- estimate$mu
      W <- estimate$W
      
      # z_all <- NULL
      # for(i in 1:M){
      #   z_all <- rbind(z_all,kronecker(diag(p_covar+1),t(z[i,])))
      # }
      z_covar <- diag(M)
      
      z_test <- z_design
      xx_test <- kronecker(I_M,x_test)
      xxz_test <- xx_test%*%z_test
      score <- t(xxz_test)%*%(yy-mu)
      xx_covar <- kronecker(I_M,rep(1,N))
      xxz_covar <- xx_covar%*%z_covar
      info <- t(xxz_test)%*%W%*%xxz_test-t(xxz_test)%*%W%*%xxz_covar%*%solve(infor_covar)%*%t(xxz_covar)%*%W%*%xxz_test
      score_test_result <- t(score)%*%solve(info)%*%score
      p_value <- 1-pchisq(as.numeric(score_test_result),df=ncol(z_test))
      return(p_value)
      
    }
    yy <- as.vector(y)
    N <- nrow(x_test)
    p_covar <- ncol(x_covar)
    M <- ncol(y)
    I_M <- diag(M)
    estimate <- Mvpoly(delta0,y,x_covar,z_design)
    mu <- estimate$mu
    W <- estimate$W
    infor_covar <- estimate$infor
    
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
    
    
    
    
    z_test <- z_design
    xx_test <- kronecker(I_M,x_test)
    xxz_test <- xx_test%*%z_test
    score <- t(xxz_test)%*%(yy-mu)
    x_covar <- cbind(1,x_covar)
    xx_covar <- kronecker(I_M,x_covar)
    xxz_covar <- xx_covar%*%z_covar
    
    info <- t(xxz_test)%*%W%*%xxz_test-t(xxz_test)%*%W%*%xxz_covar%*%solve(infor_covar)%*%t(xxz_covar)%*%W%*%xxz_test
    score_test_result <- t(score)%*%solve(info)%*%score
    p_value <- 1-pchisq(as.numeric(score_test_result),df=ncol(z_test))
    return(p_value)
  }else{
    
    sof <- "source.so"
    dyn.load(sof)
    
    
    
    nparm <- length(delta0)
    NITER <- 100
    tol   <- 1e-4
    NCOV  <- ncol(x_covar)
    ncat  <- NCOL(z_design)-1
    DEBUG     <- 2
    ret_rc    <- as.integer(1)
    ret_delta <- as.numeric(rep(-9999, nparm))
    
    ret_p     <- as.numeric(1)
    M <- ncol(y)
    n <- nrow(y)
    
    temp <- .C("score_test", as.numeric(delta0), as.integer(nparm), as.numeric(as.vector(y)),
               as.numeric(x_covar), as.numeric(x_test), as.numeric(z_design), as.integer(n), as.integer(M), 
               as.integer(ncat), as.integer(NCOV), as.integer(NITER), as.numeric(tol), 
               as.integer(DEBUG), ret_rc=ret_rc, ret_delta=ret_delta, ret_p=ret_p)
    return(temp$ret_p)
  }
  
}




score_support <- function(pxx,x_covar,z_design,y_em){
  N <- nrow(x_covar)
  p_covar <- ncol(x_covar)
  
  ###need improvement
  M <- length(y_em)/N
  ret_p <- pxx
  sof <- "getW.so"
  dyn.load(sof)
  W <- rep(0,N*M*M)
  temp <- .C("Weighted_W", as.numeric(ret_p), as.numeric(W),as.integer(N), as.integer(M))
  W.c <- temp[[2]]
  temp_mis <- .C("Weighted_W", as.numeric(y_em), as.numeric(W),as.integer(N), as.integer(M))
  W.mis <- temp_mis[[2]]
  W <- W.c-W.mis
  x_covar <- cbind(1,x_covar)
  p_col <- ncol(x_covar)
  
  sof1 <- "getXpWXp.so"
  dyn.load(sof1)
  XpWXp <- rep(0,p_col^2*M^2)
  XpWXp_C <- .C("getXpWXp",as.numeric(W),as.numeric(x_covar),as.integer(N),as.integer(M),
                as.numeric(XpWXp),as.integer(p_col))
  XpWXp <- matrix(XpWXp_C[[5]],ncol(x_covar)*M,ncol(x_covar)*M)
  
  
  sof <- "getWXp.so"
  dyn.load(sof)
  
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



