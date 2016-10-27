#' splits a kallisto experiment by gene biotype
#' @param kexp a kallisto Experiment
#' @param biotype character a selector either 'element','SINE','LINE','LTR', etc
#' @export
kexpByType<-function(kexp,biotype=NULL){

return(kexp[grep(biotype,rowRanges(kexp)$gene_biotype),])
}
