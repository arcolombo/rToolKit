#' splits a kexp into repeat level
#' @param kexp a kallistoExperiment class object assumes pre-annotated
#' @export
#' @return returns a repeat level kallisto Experiment
findRepeats<-function(kexp){

rowRanges(kexp)[which(rowRanges(kexp)$biotype_class=="repeat")]$gene_id<-rowRanges(kexp)[which(rowRanges(kexp)$biotype_class=="repeat")]$tx_id
kexp<-kexp[which(rowRanges(kexp)$biotype_class=="repeat"),]
return(kexp[!duplicated(rownames(kexp)),])

}
