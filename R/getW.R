
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
