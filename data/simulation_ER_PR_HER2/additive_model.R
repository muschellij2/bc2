#This is an example with 1000 samples;
#There are 3 covaraites in total;p.covar means the number of covariates;
#Three second stage categories;

# rm(list=ls())
# commandarg <- commandArgs(trailingOnly = T)
# i1 <- as.numeric(commandarg[1])
# print(i1)
#
# set.seed(i1)
# setwd("/home/zhangh20/breast_cancer")
logit_inver <- function(x){
  return(exp(x)/(1+exp(x)))
}

p_value_function <- function(z){
  result <- NULL
  for(i in 1:length(z)){

    result <- c(result,2*pnorm(-abs(z[i])))

  }
  return(result)
}


# five_test = function(score_value,infor){
#   globe_score <- score_value[1:4]%*%solve(infor[1:4,1:4])%*%score_value[1:4]
#   globe_p <- pchisq(as.numeric(globe_score),df=4,lower.tail = F)
#
#
#   heter_score <- score_value[2:4]%*%solve(infor[2:4,2:4])%*%score_value[2:4]
#   heter_p <- 1-pchisq(as.numeric(heter_score),df=3,lower.tail=F)
#   infor_inv = solve(infor)
#   heter_indi <- rep(0,3)
#   if(diag(infor_inv)[2]>0){
#     heter_indi[1] <-  score_value[2]/sqrt(diag(infor)[2])
#   }else{
#     heter_indi[1] <- 0
#   }
#   if(diag(infor_inv)[3]>0){
#     heter_indi[2] <-  score_value[3]/sqrt(diag(infor)[3])
#   }else{
#     heter_indi[3] <- 0
#   }
#   if(diag(infor_inv)[4]>0){
#     heter_indi[3] <-  score_value[4]/sqrt(diag(infor)[4])
#   }else{
#     heter_indi[3] <- 0
#   }
#   #heter_indi <- score_value[2:4]*sqrt(diag(infor_inv)[2:4])
#   p_value <- p_value_function(heter_indi)
#   result = c(globe_p,heter_p,p_value)
#   return(result)
# }
#




#source("mvpoly.R")








transform <- function(y,z_stanard){
  K = ncol(z_standard)
  idx <- which(y==1)
  if(length(idx)==0){
    return(rep(NA,K))
  }else{
    return(z_standard[idx,])
  }
}


error_f <- function(delta_old,delta_new){
  max(abs(delta_new-delta_old))/(max(abs(delta_old))+0.1)
}

prob_fitting <- function(delta,y_pheno,x_all,z_standard,z_all){
  result <- matrix(0,nrow=nrow(y_pheno),ncol = nrow(z_standard))
  beta <- matrix(z_all%*%delta,ncol = nrow(z_standard))
  #beta <- matrix(z_all%*%delta,nrow = ncol(z_design))
  #beta <- matrix(z_all%*%delta,nrow = nrow(z_standard))
  for(i in 1:nrow(y_pheno)){
    if(y_pheno[i,1]==1) {
      idx <- which(y_pheno[i,2:ncol(z_design)]!=888)
      if(length(idx)==0) {
        result[i,] <- exp(x_all[i,]%*%beta)/(1+sum(exp(x_all[i,]%*%beta)))
      }else{
        jdx <- apply(z_standard,1,function(t){all(t[idx]==y_pheno[i,idx+1])})
        jdx <- which(jdx==T)
        temp <- exp(x_all[i,]%*%beta[,jdx])
        result[i,jdx] <- temp/sum(temp)
      }
    }
  }

  return(result)
}


# ######complete cases ward test
 a <- c(0,1)
 b <- c(0,1)
 c <- c(0,1)
 d <- c(1:3)

 z <- as.matrix(expand.grid(a,b,c,d)) # orig
 #z <- as.matrix(expand.grid(a,b))
 z_standard <- z
#
 y_standard <- diag(nrow(z_standard))
#
 #this z_design matrix is the second stage matrix
 z_design <- cbind(1,z)
 M <- nrow(z_design)
 p.covar = 2
 z_all <- NULL
 z_all_temp <- NULL
 for(i in 1:M){
   z_all_temp <- rbind(z_all_temp,kronecker(diag(p.covar),t(z_design[i,])))
 }
#
 z_all <- matrix(0,nrow = M*(p.covar+1),ncol= M+p.covar*ncol(z_design))
 for(i in 1:M){
   z_all[1+(i-1)*(p.covar+1),i] <- 1
 }
 for(i in 1:M){
   z_all[(2+(i-1)*(p.covar+1)):(i*(p.covar+1)),(M+1):ncol(z_all)] <-
     z_all_temp[(1+(i-1)*p.covar):(i*p.covar),]
 }
 # for(i in 1:(M)){
 #   temp <- rep(0,ncol(z_all))
 #   temp[i] <- 1
 #   z_all[1+(i-1)*(p.covar+1),] = temp
 # }
 K <- ncol(z_design)

 # z <- kronecker(diag(2),z)
 theta_intercept <- rep(1,M)
 theta_test <- c(0,rep(0,K-1))
 theta_covar <- rep(rep(0.5,K),p.covar-1)

 #this theta is the true value
 theta <- c(theta_intercept,theta_test,theta_covar)

 #this is the true beta
 beta <- z_all%*%theta
 beta <- matrix(beta,nrow=p.covar+1)
 simulationtimes <- 1
 #delta_result <- matrix(0,nrow=simulationtimes,ncol = length(theta))
 #infor_result <- matrix(0,nrow=simulationtimes*ncol(z_all),ncol = ncol(z_all))

 result_wald_complete = matrix(0,nrow=simulationtimes,ncol=(ncol(z_design)+ncol(z_design)^2))
 p_wald_complete = matrix(0,nrow=simulationtimes,5)




#for(simulation in 1:simulationtimes){
  #print(paste0("we are in",simulation,"th simulation"))

  #alpha <- c(0,rep(1,length(beta)-1))

  n <- 10000
  x <-  matrix(rnorm(p.covar*n),nrow = n)
  x_test <- x[,1]
  x_covar <- x[,2:ncol(x)]
  x_all <- cbind(1,x) ##adding the intercept into the model


  predictor <- x_all%*%(beta)
  predictor <- cbind(0,predictor)
  #predictor <- sweep(predictor,2,alpha,"+")
  p <- exp(predictor)
  sp <- rowSums(p)
  #this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
  p <- sweep(p,1,sp,"/")
  y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))

  y <- y[,-1] # this is the y matrix for the model


  ###creating the y phenotype file which user will use
  y_pheno <- matrix(0,nrow = nrow(y), ncol = K)
  y_control <- rep(0,ncol(y))
  idx <- apply(y,1,function(x){all(x==y_control)})
  idx <- !idx

  y_pheno[idx,1] <- 1
  y_pheno[,2:K] <- t(apply(y,1,function(x){transform(x,z_standard)}))
  idx <- which(idx == T)
  ###generating missing data
   missing_rate <- 0.2
   for(i in 2:K){
     eval(parse(text = paste0("y_pheno[sample(idx,missing_rate*length(idx)),",i,"]=888")))
   }



  delta0 <- c(theta_intercept,theta_test,theta_covar)



  epi <- 1e-04
  #delta_old <- rep(0,length(delta0))
  delta_old <- theta
  ##EM algorithm
  for(niter in 1:100){

    y_em <- prob_fitting(delta_old,y_pheno,x_all,z_standard,z_all)
    result <-  Mvpoly(delta_old,y_em,x,z_design)
    delta_new <-result[[1]]
    error_f(delta_old,delta_new)
    if(error_f(delta_old,delta_new)<epi){

      print("Converge!!!!!")
      break
    }
    delta_old <- delta_new
    print(niter)

  }
  delta_temp = delta_new[ (length(theta_intercept)+1):
                            + (length(theta_intercept)+ncol(z_design))]

  infor_obs_inv = solve(result[[2]])[(length(theta_intercept)+1):
                                       + (length(theta_intercept)+ncol(z_design)),(length(theta_intercept)+1):
                                       + (length(theta_intercept)+ncol(z_design))]
  temp_result = c(delta_temp,as.vector(infor_obs_inv))

  result_wald_complete[simulation,] <- temp_result
  p_wald_complete = five_test(delta_temp,infor_obs_inv)

#}


data = list(y_pheno=y_pheno,x=x)
save(data,file = "./data/simulation_ER_PR_HER2/simulated_data.Rdata")








###score test complete cases
a <- c(0,1)
b <- c(0,1)
c <- c(0,1)
p.covar <- 3
z <- as.matrix(expand.grid(a,b,c)) # orig
#z <- as.matrix(expand.grid(a,b))
z_standard <- z

y_standard <- diag(nrow(z_standard))

#this z_design matrix is the second stage matrix
z_design <- cbind(1,z)
M <- nrow(z_design)
p.covar <- 2
z_all <- NULL
z_all_temp <- NULL
for(i in 1:M){
  z_all_temp <- rbind(z_all_temp,kronecker(diag(p.covar),t(z_design[i,])))
}

z_all <- matrix(0,nrow = M*(p.covar+1),ncol= M+p.covar*ncol(z_design))
for(i in 1:M){
  z_all[1+(i-1)*(p.covar+1),i] <- 1
}
for(i in 1:M){
  z_all[(2+(i-1)*(p.covar+1)):(i*(p.covar+1)),(M+1):ncol(z_all)] <-
    z_all_temp[(1+(i-1)*p.covar):(i*p.covar),]
}
# for(i in 1:(M)){
#   temp <- rep(0,ncol(z_all))
#   temp[i] <- 1
#   z_all[1+(i-1)*(p.covar+1),] = temp
# }
K <- ncol(z_design)

# z <- kronecker(diag(2),z)
theta_intercept <- rep(1,M)
theta_test <- c(0,rep(0,K-1))
theta_covar <- rep(rep(0.5,K),p.covar)

#this theta is the true value
theta <- c(theta_intercept,theta_covar)

#this is the true beta
beta <- z_all%*%theta
beta <- matrix(beta,nrow=(p.covar+1))
simulationtimes <- 1
delta_result <- matrix(0,nrow=simulationtimes,ncol = length(theta))
infor_result <- matrix(0,nrow=simulationtimes*ncol(z_all),ncol = ncol(z_all))


n <- 50000
x <-  matrix(rnorm(p.covar*n),nrow = n)

#x_test <- x_test_all[,simulation]
x_covar <- x
x_all <- cbind(1,x_covar) ##adding the intercept into the model


predictor <- x_all%*%(beta)
predictor <- cbind(0,predictor)
#predictor <- sweep(predictor,2,alpha,"+")
p <- exp(predictor)
sp <- rowSums(p)
#this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
p <- sweep(p,1,sp,"/")
y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))

y <- y[,-1] # this is the y matrix for the model


###creating the y phenotype file which user will use
y_pheno <- matrix(0,nrow = nrow(y), ncol = 4)
y_control <- rep(0,ncol(y))
idx <- apply(y,1,function(x){all(x==y_control)})
idx <- !idx

y_pheno[idx,1] <- 1
y_pheno[,2:4] <- t(apply(y,1,function(x){transform(x,z_standard)}))
idx <- which(idx == T)
###generating missing data
#missing_rate <- 0.2
#y_pheno[sample(idx,missing_rate*length(idx)),2] =888
#y_pheno[sample(idx,missing_rate*length(idx)),3] =888
#y_pheno[sample(idx,missing_rate*length(idx)),4] =888




delta0 <- c(theta_intercept,theta_covar)

epi <- 1e-04
#delta_old <- rep(0,length(delta0))
delta_old <- delta0
##EM algorithm
for(niter in 1:100){

  y_em <- prob_fitting(delta_old,y_pheno,x_all,z_standard,z_all)
  result <-  Mvpoly(delta_old,y_em,x,z_design)
  delta_new <-result[[1]]
  error_f(delta_old,delta_new)
  if(error_f(delta_old,delta_new)<epi){

    print("Converge!!!!!")
    break
  }
  delta_old <- delta_new
  print(niter)

}
pxx <- result[[3]]
score_support_result <- score_support(pxx,x_covar,z_design,y_em)

simulationtimes <- 1
score_result <- matrix(0,nrow = simulationtimes,4)
infor_result <- matrix(0,nrow=simulationtimes*4,4)
p_result_global <- rep(0,simulationtimes)
p_result_heter <- rep(0,simulationtimes)
p_result_ind <- matrix(0,nrow=simulationtimes,3)
x_test_all <- matrix(rnorm(n*simulationtimes),nrow=n)
p_score_complete = matrix(0,simulationtimes,5)
result_score_complete = matrix(0,nrow=simulationtimes,ncol=(ncol(z_design)+ncol(z_design)^2))


for(simulation in 1:simulationtimes){
  print(simulation)
  snpvalue <- x_test_all[,simulation]
  temp_score_result  <- score_test_mis(y_em,snpvalue,score_support_result)
  score_value <- temp_score_result[[1]]
  score_result[simulation,] <-  score_value
  infor <- temp_score_result[[2]]
  infor_result[4*(simulation-1)+(1:4),] <- infor
  result_score_complete[simulation,] = c(score_value,as.vector(infor))
  p_score_complete[simulation,] = five_test(score_value,infor)

}










#####
######incomplete cases ward test
a <- c(0,1)
b <- c(0,1)
c <- c(0,1)
p.covar <- 3
z <- as.matrix(expand.grid(a,b,c)) # orig
#z <- as.matrix(expand.grid(a,b))
z_standard <- z

y_standard <- diag(nrow(z_standard))

#this z_design matrix is the second stage matrix
z_design <- cbind(1,z)
M <- nrow(z_design)
z_all <- NULL
z_all_temp <- NULL
for(i in 1:M){
  z_all_temp <- rbind(z_all_temp,kronecker(diag(p.covar),t(z_design[i,])))
}

z_all <- matrix(0,nrow = M*(p.covar+1),ncol= M+p.covar*ncol(z_design))
for(i in 1:M){
  z_all[1+(i-1)*(p.covar+1),i] <- 1
}
for(i in 1:M){
  z_all[(2+(i-1)*(p.covar+1)):(i*(p.covar+1)),(M+1):ncol(z_all)] <-
    z_all_temp[(1+(i-1)*p.covar):(i*p.covar),]
}
# for(i in 1:(M)){
#   temp <- rep(0,ncol(z_all))
#   temp[i] <- 1
#   z_all[1+(i-1)*(p.covar+1),] = temp
# }
K <- ncol(z_design)

# z <- kronecker(diag(2),z)
theta_intercept <- rep(1,M)
theta_test <- c(0,rep(0,K-1))
theta_covar <- rep(rep(0.5,K),p.covar-1)

#this theta is the true value
theta <- c(theta_intercept,theta_test,theta_covar)

#this is the true beta
beta <- z_all%*%theta
beta <- matrix(beta,nrow=p.covar+1)
simulationtimes <- 1
#delta_result <- matrix(0,nrow=simulationtimes,ncol = length(theta))
#infor_result <- matrix(0,nrow=simulationtimes*ncol(z_all),ncol = ncol(z_all))

result_wald_incomplete = matrix(0,nrow=simulationtimes,ncol=(ncol(z_design)+ncol(z_design)^2))
p_wald_incomplete = matrix(0,nrow=simulationtimes,5)

for(simulation in 1:simulationtimes){
  print(paste0("we are in",simulation,"th simulation"))

  #alpha <- c(0,rep(1,length(beta)-1))

  n <- 50000
  x <-  matrix(rnorm(p.covar*n),nrow = n)
  x_test <- x[,1]
  x_covar <- x[,2:ncol(x)]
  x_all <- cbind(1,x) ##adding the intercept into the model


  predictor <- x_all%*%(beta)
  predictor <- cbind(0,predictor)
  #predictor <- sweep(predictor,2,alpha,"+")
  p <- exp(predictor)
  sp <- rowSums(p)
  #this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
  p <- sweep(p,1,sp,"/")
  y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))

  y <- y[,-1] # this is the y matrix for the model


  ###creating the y phenotype file which user will use
  y_pheno <- matrix(0,nrow = nrow(y), ncol = 4)
  y_control <- rep(0,ncol(y))
  idx <- apply(y,1,function(x){all(x==y_control)})
  idx <- !idx

  y_pheno[idx,1] <- 1
  y_pheno[,2:4] <- t(apply(y,1,function(x){transform(x,z_standard)}))
  idx <- which(idx == T)
  ###generating missing data
  missing_rate <- 0.2
  y_pheno[sample(idx,missing_rate*length(idx)),2] =888
  y_pheno[sample(idx,missing_rate*length(idx)),3] =888
  y_pheno[sample(idx,missing_rate*length(idx)),4] =888




  delta0 <- c(theta_intercept,theta_test,theta_covar)



  epi <- 1e-04
  #delta_old <- rep(0,length(delta0))
  delta_old <- theta
  ##EM algorithm
  for(niter in 1:100){

    y_em <- prob_fitting(delta_old,y_pheno,x_all,z_standard,z_all)
    result <-  Mvpoly(delta_old,y_em,x,z_design)
    delta_new <-result[[1]]
    error_f(delta_old,delta_new)
    if(error_f(delta_old,delta_new)<epi){

      print("Converge!!!!!")
      break
    }
    delta_old <- delta_new
    print(niter)

  }
  infor_mis_c <- infor_mis(y_em,x,z_design)
  infor_obs <- result[[2]]-infor_mis_c
  delta_temp = delta_new[ (length(theta_intercept)+1):
                            (length(theta_intercept)+ncol(z_design))]
  infor_obs_inv = solve(infor_obs)[(length(theta_intercept)+1):
                                     (length(theta_intercept)+ncol(z_design)),(length(theta_intercept)+1):
                                     (length(theta_intercept)+ncol(z_design))]
  temp_result = c(delta_temp,as.vector(infor_obs_inv))

  result_wald_incomplete[simulation,] <- temp_result
  p_wald_incomplete = five_test(delta_temp,infor_obs_inv)

}







###score test incomplete cases
a <- c(0,1)
b <- c(0,1)
c <- c(0,1)
d <- c(1:3)
K <- 4
z <- as.matrix(expand.grid(a,b,c,d)) # orig
#z <- as.matrix(expand.grid(a,b))
z_standard <- z

y_standard <- diag(nrow(z_standard))

#this z_design matrix is the second stage matrix
z_design <- cbind(1,z)
M <- nrow(z_design)
p.covar <- 2
z_all <- NULL
z_all_temp <- NULL
for(i in 1:M){
  z_all_temp <- rbind(z_all_temp,kronecker(diag(p.covar),t(z_design[i,])))
}

z_all <- matrix(0,nrow = M*(p.covar+1),ncol= M+p.covar*ncol(z_design))
for(i in 1:M){
  z_all[1+(i-1)*(p.covar+1),i] <- 1
}
for(i in 1:M){
  z_all[(2+(i-1)*(p.covar+1)):(i*(p.covar+1)),(M+1):ncol(z_all)] <-
    z_all_temp[(1+(i-1)*p.covar):(i*p.covar),]
}
# for(i in 1:(M)){
#   temp <- rep(0,ncol(z_all))
#   temp[i] <- 1
#   z_all[1+(i-1)*(p.covar+1),] = temp
# }
K <- ncol(z_design)

# z <- kronecker(diag(2),z)
theta_intercept <- rep(1,M)
theta_test <- c(0,rep(0,K-1))
theta_covar <- rep(rep(0.5,K),p.covar)

#this theta is the true value
theta <- c(theta_intercept,theta_covar)

#this is the true beta
beta <- z_all%*%theta
beta <- matrix(beta,nrow=(p.covar+1))
simulationtimes <- 1
delta_result <- matrix(0,nrow=simulationtimes,ncol = length(theta))
infor_result <- matrix(0,nrow=simulationtimes*ncol(z_all),ncol = ncol(z_all))


n <- 50000
x <-  matrix(rnorm(p.covar*n),nrow = n)

#x_test <- x_test_all[,simulation]
x_covar <- x
x_all <- cbind(1,x_covar) ##adding the intercept into the model


predictor <- x_all%*%(beta)
predictor <- cbind(0,predictor)
#predictor <- sweep(predictor,2,alpha,"+")
p <- exp(predictor)
sp <- rowSums(p)
#this standarize the probability sum into 1,since logit model:pr(D=1|predictor)=exp(predictor)/(1+exp(predictor))
p <- sweep(p,1,sp,"/")
y <- t(apply(p,1,function(x){rmultinom(1,1,x)}))

y <- y[,-1] # this is the y matrix for the model


###creating the y phenotype file which user will use
y_pheno <- matrix(0,nrow = nrow(y), ncol = 4)
y_control <- rep(0,ncol(y))
idx <- apply(y,1,function(x){all(x==y_control)})
idx <- !idx

y_pheno[idx,1] <- 1
y_pheno[,2:(K+1)] <- t(apply(y,1,function(x){transform(x,z_standard)}))
idx <- which(idx == T)
###generating missing data
missing_rate <- 0.2
y_pheno[sample(idx,missing_rate*length(idx)),2] =888
y_pheno[sample(idx,missing_rate*length(idx)),3] =888
y_pheno[sample(idx,missing_rate*length(idx)),4] =888
y_pheno[sample(idx,missing_rate*length(idx)),5] = 888




delta0 <- c(theta_intercept,theta_covar)

epi <- 1e-04
#delta_old <- rep(0,length(delta0))
delta_old <- delta0
##EM algorithm
for(niter in 1:100){

  y_em <- prob_fitting(delta_old,y_pheno,x_all,z_standard,z_all)
  result <-  Mvpoly(delta_old,y_em,x,z_design)
  delta_new <-result[[1]]
  error_f(delta_old,delta_new)
  if(error_f(delta_old,delta_new)<epi){

    print("Converge!!!!!")
    break
  }
  delta_old <- delta_new
  print(niter)

}
pxx <- result[[3]]
score_support_result <- score_support(pxx,x_covar,z_design,y_em)

simulationtimes <- 1
score_result <- matrix(0,nrow = simulationtimes,4)
infor_result <- matrix(0,nrow=simulationtimes*4,4)
p_result_global <- rep(0,simulationtimes)
p_result_heter <- rep(0,simulationtimes)
p_result_ind <- matrix(0,nrow=simulationtimes,3)
x_test_all <- matrix(rnorm(n*simulationtimes),nrow=n)
p_score_incomplete = matrix(0,simulationtimes,5)
result_score_incomplete = matrix(0,nrow=simulationtimes,ncol=(ncol(z_design)+ncol(z_design)^2))


for(simulation in 1:simulationtimes){
  print(simulation)
  snpvalue <- x_test_all[,simulation]
  temp_score_result  <- score_test_mis(y_em,snpvalue,score_support_result)
  score_value <- temp_score_result[[1]]
  score_result[simulation,] <-  score_value
  infor <- temp_score_result[[2]]
  infor_result[4*(simulation-1)+(1:4),] <- infor
  result_score_incomplete[simulation,] = c(score_value,as.vector(infor))
  p_score_incomplete[simulation,] = five_test(score_value,infor)

}



final_result = list(p_wald_complete,result_wald_complete,p_score_complete,result_score_complete,
                    p_wald_incomplete,result_wald_incomplete,p_score_incomplete,result_score_incomplete)


save(final_result,file = paste0("/data/zhangh20/breast_cancer/simulation1/final_result",i1))



