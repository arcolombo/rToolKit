#' @title finds the non monotonic nominal ordering group
#' @description for a given stage the ordering is selected
#' @param mt a matrix of collapse bundles of repeats
#' @param globalMax the selected global max column number
#' @param globalMin the selected global min column number
#' @param globalMid the column number of desired mid
#' @return a matrix group that fits into this nominal group ordering
#' @export
orderMonotonicity<-function(mt,globalMax=stageId1,globalMin=stageId2,globalMid=stageId3){
  
  id1<- which(mt[,globalMin]<mt[,globalMid])
  tt<-mt[id1,]
  
  stopifnot(all(tt[,globalMin]<tt[,globalMid])==TRUE)
  id2<-which(tt[,globalMid]<tt[,globalMax])
  t2<-tt[id2,]
  if(is.null(dim(t2))==TRUE){
  t2<-t(as.data.frame(t2))
  }
  stopifnot(all(t2[,globalMax]>t2[,globalMid]))
  
  stopifnot(all(t2[,globalMid]>t2[,globalMin])==TRUE)
  stopifnot(all(t2[,globalMax]>t2[,globalMid])==TRUE)
  stopifnot(all(t2[,globalMax]>t2[,globalMin])==TRUE)
  return(t2)

} #{{{ main
