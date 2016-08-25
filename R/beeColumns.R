#' @title beeswarm plot of expression columns
#' @description plots beeswarm of all repeats by clonal stage with options for normalized CPM , note that TPM is not normalized.
#' @param kexp a kexp
#' @param patientID a character to select patient info
#' @param repeat_biotype a character either 'element','LTR','SINE', etc
#' @param selected how many byMads to select
#' @param what  is count or tpm
#' @param how   norm method TMM CQN or none (raw est_counts)
#' @import beeswarm
#' @import arkas
#' @export
beeColumns<-function(kexp,patientID=NULL,repeat_biotype=NULL,selected=20,what=c("tpm","counts"), how=c("TMM","CQN","none")){
  
  
  rep.names<-rownames(byMad(counts(kexp),k=selected))
  phenodData<-as.data.frame(pData(kexp))
  what <- match.arg(what,c("counts","tpm"))
  how<-match.arg(how,c("TMM","CQN","none"))
  if(what=="counts"){
      if(how=="none"){
  dataMat<-as.data.frame(counts(kexp)[rep.names,])
   } else if(how=="TMM"){
  stopifnot(is.null(assays(kexp)$cpm_library_normalized)==FALSE)
  dataMat<-as.data.frame(assays(kexp)$cpm_library_normalized[rep.names,])
   } else  if(how=="CQN") {
   #FIX ME:  add CQN data 
   message("CQN is not supported yet..")
  # stopifnot(is.null(assays(kexp)....==FALSE)
  # datamat<-...
   }
  } else {
  dataMat<-as.data.frame(tpm(kexp)[rep.names,])
  }
  beeswarm(asinh(dataMat), pch=16,xlab=paste0(patientID," ",repeat_biotype),
  ylab=paste0(toupper(what), " Norm. ",how))
 
}
