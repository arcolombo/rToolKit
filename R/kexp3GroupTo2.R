#' splits kexp by group  pHSC , LSC, Blast
#' @param kexp a kallisto experiment
#' @param comparison  pHSC LSC etc
#' @param control  not comparions
#' @export 
kexp3GroupTo2 <-function(kexp,comparison=NULL,newComparison=NULL){
  id1<-grep(comparison,colnames(kexp))
  colnames(kexp)[grep(comparison,colnames(kexp))]<-paste0(colnames(kexp)[grep(comparison,colnames(kexp))],"_1")
 
  colnames(kexp)<-gsub(comparison,newComparison,colnames(kexp))
  pData(kexp)$ID<-colnames(kexp)
 
  group = sapply(strsplit(pData(kexp)$ID,"_"),function(x) x[1])
  design<-model.matrix(~group)
  rownames(design)<-pData(kexp)$ID
  design[grepl(comparison,rownames(design)),2]<-1


  metadata(kexp)$design<-design
  print(metadata(kexp)$design)
return(kexp)
}
