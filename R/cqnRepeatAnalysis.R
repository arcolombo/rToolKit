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
cqnRepeatAnalysis<-function(fullKexp,cqnRepeat=NULL,inversion=TRUE,comparison=comparison,control=control){
  if(is.null(cqnRepeat)==TRUE){
  cqn.full<-cqnDE(fullKexp,comparison=comparison,control=control)
  } 
  if(inversion==TRUE){
  Repeat.counts<-cqn.full[["RPKM.cqn.log2"]]
  Repeat.counts<-Repeat.counts^2
  } 
  Repeat.counts<-Repeat.counts[!grepl("^ERCC",rownames(Repeat.counts)),]
  Repeat.counts<-Repeat.counts[!grepl("^ENST",rownames(Repeat.counts)),]
  fullTags<-rownames(cqnRepeat["topTags"][[1]])
  ##now to compar
  rep<-findRepeats(fullKexp)
  
  onlyRepeats<-cqnDE(rep,comparison=comparison,control=control)
  subsetTags<-rownames(onlyRepeats["topTags"][[1]])
  tester<-kexp2Group(fullKexp,comparison=comparison,control=control)
  Heatmap(asinh(counts(tester)[rownames(counts(tester))%in%subsetTags,]),main=paste0("CQN Normalization Only Using CPM Repeats") )
  Heatmap(asinh(tpm(tester)[rownames(counts(tester))%in%subsetTags,]),main="CQN Normalization TPM Only Repeats"  )
 
  ###look at CQN norm top genes with limma and compare with RUV should show same trend



}


