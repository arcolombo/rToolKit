#' @title cqnDE analysis on repeats, major differences in performance when CQN is called on transcripts and/or called on repeat elements. ALUs have higher GC conetnet and may have bias
#' @description this methods is downstream of cqnDE where it looks at repeat subset expression from the cqnDE entire transcript set.  we are not sure if the cqn method has a bias with regard to repeats, and we are not sure if repeats with sequence repetitive nature can bias the cqn method when **only** normalizing with respect to repeats.   so cqnDE is called with the entire transcript set, and the repeats are merely subset out; contrast to calling cqn on repeats only and analyzing the result.
#' @param fullKexp a full kallistoExperiment, by Stage is ok
#' @param cqnRepeat  a data matrix of cqnDE result initially from a full expression log2 scale of all transcripts, a full normalization
#' @param inversion  boolean , if true then the log2 inverse is performed
#' @import ComplexHeatmap
#' @importFrom quantreg rq
#' @import cqn
#' @export
#' @return  returns some images on comparing repeat CQN method
cqnRepeatAnalysis<-function(fullKexp,cqnRepeat=NULL,inversion=TRUE,comparison=comparison,control=control,whichDelta=c("delta1","delta2")){
  if(is.null(cqnRepeat)==TRUE){
  cqn.full<-cqnDE(fullKexp,comparison=comparison,control=control)
  } 
  whichDelta=match.arg(whichDelta,c("delta1","delta2"))
  if(inversion==TRUE){
  Repeat.counts<-cqn.full[["RPKM.cqn.log2"]]
  Repeat.counts<-Repeat.counts^2
  } 
  Repeat.counts<-Repeat.counts[!grepl("^ERCC",rownames(Repeat.counts)),]
  Repeat.counts<-Repeat.counts[!grepl("^ENST",rownames(Repeat.counts)),]
  fullTags<-rownames(cqn.full["topTags"][[1]])
  fullDE<-as.data.frame(cqn.full["topTags"][1])
  colnames(fullDE)<-c("length","gccontent","logFC","logCPM","LR","PValue","FDR")
  if(whichDelta=="delta1"){
  plotFrequency(fullKexp,topNames=fullTags,topDE=fullDE,whichDelta="delta1")
  beeswarm(fullDE$logFC, main=expression(paste(Delta,"(pHSC,LSC) CQN")),xlab=paste0("LSC-pHSC p.val 0.05"))
  par(new=TRUE)
  boxplot(fullDE$logFC,medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  }
  if(whichDelta=="delta2"){
  plotFrequency(fullKexp,topNames=fullTags,topDE=fullDE,whichDelta="delta2")
  beeswarm(fullDE$logFC, main=expression(paste(Delta,"(LSC,Blast) CQN")),xlab=paste0("Blast-LSC p.val 0.05"))
  par(new=TRUE)
  boxplot(fullDE$logFC,medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  }
##heatmapping the selected GC normalized DE repeat list in terms of the original kexp counts follows the pattern
###yyz<-(counts(fullKexp))[rownames(counts(fullKexp))%in%fullTags,]
###Heatmap(log(1+yyz))

##heatmapping the selected DE repeats in terms of the normalized CQN counts shows an inverse pattern with LSC clearly higher than pHSC in raw counts.  however the overall DE and logFC follows the pattern, but the LSC expression pattern inverses
##zyy<-Repeat.counts[rownames(Repeat.counts)%in%fullTags,]
###Heatmap(log(1+zzy))
##one reason why LSC are upregulated in the cqn counts is because it is determined by gc content and length adjusted values.
###


}  ###main


