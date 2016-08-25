#' beeswarm patient info plotting cpm or tpm of patient specific, repeat specififc of three interesting clonal processes.
#' @param kexp a kexp
#' @param patientID a character to select patient info
#' @param repeat_biotype a character either 'element','LTR','SINE', etc
#' @param how  norm type TMM CQN none 
#' @import beeswarm
#' @import arkas
#' @export
beePatient<-function(kexp, patientID=NULL, repeat_biotype=NULL,selected=20,how=c("TMM","CQN","none")){
 ##FIX ME: add plotting for selected!=20  
  ###need to select top k MADs or SDs for the selection
  how<-match.arg(how,c("TMM","CQN","none"))
 if(how=="none"){
 rep.names<-rownames(byMad(counts(kexp),k=selected))
 } else if(how=="TMM"){
  rep.names<-rownames(byMad(assays(kexp)$cpm_library_normalized,k=selected))
 } else if(how=="CQN"){
  message("CQN not yet..")
 }
 
 phenoData<-as.data.frame(pData(kexp))

#rpts is built from collapseTranscripts and filtering out only repeats
# plot 'em: counts
  if(how=="none"){
 op<- par(mfrow=c(5,4))
  for ( name in rep.names){
   getBeeswarm(name, dataMat=counts(kexp), what="count", phenoData,patientID,repeat_biotype,how)
   }
  par(op)
  readkey()
 op2<- par(mfrow=c(5,4))
  for ( name in rep.names) {
  getBeeswarm(name, dataMat=counts(kexp), what="tpm", phenoData,patientID,repeat_biotype, how)
  }
 par(op2)
 readkey()
 } else if(how=="TMM"){
 op<- par(mfrow=c(5,4))
  for ( name in rep.names){
   getBeeswarm(name, dataMat=assays(kexp)$cpm_library_normalized, what="count", phenoData,patientID,repeat_biotype,how)
   }
  par(op)
  readkey()
 op2<- par(mfrow=c(5,4))
  for ( name in rep.names) {
  getBeeswarm(name, dataMat=counts(kexp), what="tpm", phenoData,patientID,repeat_biotype, how)
  }
 par(op2)
 readkey()
  }

} #main
