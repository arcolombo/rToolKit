#' orders selection by sd
#' @param x a matrix like object
#' @param k the number of top selected sd
#' @export
#' @importFrom stats sd
bySd <- function(x, k=500) { 
   sds<-vector() 
  for(i in 1:nrow(x)){
   sds[i]<-sd(x[i,],na.rm=TRUE)
   }
  x[rev(order(sds))[seq_len(k)],]
}

