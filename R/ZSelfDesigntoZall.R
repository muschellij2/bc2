ZSelfDesigntoZall <- function(z.design,x){
  p.col <- ncol(x)
  M <- nrow(z.design)
  z.all <- NULL
  z.all.temp <- NULL
  for(i in 1:M){
    z.all.temp <- rbind(z.all.temp,kronecker(diag(p.col),t(z.design[i,])))
  }

  z.all <- matrix(0,nrow = M*(p.col+1),ncol= M+p.col*ncol(z.design))
  for(i in 1:M){
    z.all[1+(i-1)*(p.col+1),i] <- 1
  }
  for(i in 1:M){
    z.all[(2+(i-1)*(p.col+1)):(i*(p.col+1)),(M+1):ncol(z.all)] <-
      z.all.temp[(1+(i-1)*p.col):(i*p.col),]
  }

}
