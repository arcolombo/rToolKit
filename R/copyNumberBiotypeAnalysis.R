#' @title analyzes the biotypes of the copy number elements
#' @description using the copy number information, analyzes the biotypes of each copy number and ranks them
#' @param kexp a kallistoExperiment 
#' @import TxDbLite
#' @export
#' @return null 
copyNumberBiotypeAnalysis<-function(kexp){
require(RepDbLite.Hsapiens.2103)
 repTx<-transcripts(RepDbLite("RepDbLite.Hsapiens.2103"))
 
 data("repeatCopyNumber.hg19",package="TxDbLite")
cn<-repeatCopyNumber.hg19
total<-sum(cn$count)
 repDF<-as.data.frame(repTx) 
 repDF<-repDF[,c(8,13,14)]
 ff<-split(repDF,repDF$gene_biotype)
 tt<-lapply(ff,function(x) sum(x$copyNumber))
 gene_perc<-lapply(tt,function(x) (x/total)*100)

 gg<-split(repDF,repDF$tx_biotype)
 mm<-lapply(gg,function(x) sum(x$copyNumber))
 tx_perc<-lapply(mm,function(x) (x/total)*100)
 print(tx_perc)
 print(gene_perc)
 return(NULL)
}
