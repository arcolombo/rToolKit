#' @title adds repeat annotations to gene level analysis
#' @description allows for gene level analsysis competitive against genes and not only transcripts
#' @import arkas
#' @export
#' @return a kexp
fixRepeatRanges<-function(kexp){

 rowRanges(kexp)[which(rowRanges(kexp)$biotype_class=="repeat")]$gene_id<-rowRanges(kexp)[which(rowRanges(kexp)$biotype_class=="repeat")]$tx_id
 return(kexp)
}
