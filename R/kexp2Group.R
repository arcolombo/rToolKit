#' splits kexp by group  pHSC , LSC, Blast
#' @param kexp a kallisto experiment
#' @param comparison  pHSC LSC etc
#' @param control  not comparions
#' @export 
kexp2Group <-function(kexp,comparison=NULL,control=NULL){
  id1<-grep(comparison,colnames(kexp))
  id2<-grep(control,colnames(kexp))
  kexp<-kexp[,c(id1,id2)] 
  group = sapply(strsplit(pData(kexp)$ID,"_"),function(x) x[1])
  design<-model.matrix(~group)
  rownames(design)<-pData(kexp)$ID
  design[grepl(comparison,rownames(design)),2]<-1
  design[grepl(control,rownames(design)),2]<-0

  metadata(kexp)$design<-design
  print(metadata(kexp)$design)
return(kexp)
}
