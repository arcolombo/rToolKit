#' runs a beeswarm plot
#' @import beeswarm
#' @export
getBeeswarm <- function(name, dataMat, what=c("tpm", "count" ),phenoData, patientID,repeat_biotype) {

  if(is.null(repeat_biotype)==TRUE){
  repeat_biotype<-"Repetitive Elements"
  }
  what <- match.arg(what,c("counts","tpm"))
  beeswarm(dataMat[name, ] ~ phenoData$ID, col=c("green","red","blue"),
           pch=16, xlab=paste0(patientID," ",repeat_biotype),
           ylab=paste(toupper(what), "for", name))

}
