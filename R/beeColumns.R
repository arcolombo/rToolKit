#' plots beeswarm of all repeats by clonal stage
#' @param kexp a kexp
#' @param patientID a character to select patient info
#' @param repeat_biotype a character either 'element','LTR','SINE', etc
#' @param selected how many byMads to select
#' @param what  is count or tpm
#' @import beeswarm
#' @import arkas
#' @export
beeColumns<-function(kexp,patientID=NULL,repeat_biotype=NULL,selected=20,what=c("tpm","counts")){
  
  
  rep.names<-rownames(byMad(counts(kexp),k=selected))
  phenodData<-as.data.frame(pData(kexp))
  what <- match.arg(what,c("counts","tpm"))
  if(what=="counts"){
  dataMat<-as.data.frame(counts(kexp)[rep.names,])
  } else {
  dataMat<-as.data.frame(tpm(kexp)[rep.names,])
  }
  beeswarm(asinh(dataMat), pch=16,xlab=paste0(patientID," ",repeat_biotype),
  ylab=paste0(toupper(what), "for ",patientID))
 
}
