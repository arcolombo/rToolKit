#' beeswarm patient info plotting cpm or tpm of patient specific, repeat specififc of three interesting clonal processes.
#' @param kexp a kexp
#' @param patientID a character to select patient info
#' @param repeat_biotype a character either 'element','LTR','SINE', etc
#' 
#' @import beeswarm
#' @import arkas
#' @export
beePatient<-function(kexp, patientID=NULL, repeat_biotype=NULL,selected=20){
 ##FIX ME: add plotting for selected!=20  
  ###need to select top k MADs or SDs for the selection
 rep.names<-rownames(byMad(counts(kexp),k=selected))
 #stopifnot(is.null(phenoDatas)!=TRUE)
 
 phenoData<-as.data.frame(pData(kexp))

#rpts is built from collapseTranscripts and filtering out only repeats
# plot 'em: counts
 op<- par(mfrow=c(5,4))
  for ( name in rep.names){
   getBeeswarm(name, dataMat=counts(kexp), what="count", phenoData,patientID,repeat_biotype)
   }
  par(op)
  readkey()
 op2<- par(mfrow=c(5,4))
  for ( name in rep.names) {
  getBeeswarm(name, dataMat=counts(kexp), what="tpm", phenoData,patientID,repeat_biotype)
  }
 par(op2)
 readkey()
}
