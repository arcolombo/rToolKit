#' @title draws a heatmap by MAD
#' @description a general function for heatmapping data by MAD
#' @param data a matrix of data values
#' @param kexp a kallistoExperiment by any stage, or could be repeats, general use 
#' @param topAnnoFactors a factor of top annotation bars
#' @import ComplexHeatmap
#' @import stats
#' @import arkas
#' @export
#' @return returns images
heatmapByMad<-function(data,kexp,topAnnoFactors=NULL,selectK=NULL,byWhat=byWhat){
  if(is.null(selectK)==TRUE){
  selectK<-nrow(data)/2
  } else {
  selectK<-100
  }
  if(is.null(topAnnoFactors)==FALSE){
  te<-HeatmapAnnotation(as.data.frame(topAnnoFactors))
  inputMad<-byMad(data,k=selectK)
  typed<-Heatmap(inputMad,column_title=paste0("byMad Top ",selectK," ",byWhat), top_annotation=te)
  } else {
  typed<-Heatmap(inputMad,column_title=paste0("byMad Top ",selectK," ",byWhat))
  }
  typed2<-Heatmap(rowRanges(kexp)[rownames(inputMad)]$tx_biotype,name="tx_biotype",width=unit(5,"mm")) 
  typed3<-Heatmap(rowRanges(kexp)[rownames(inputMad)]$gene_biotype, name="gene_biotype",width=unit(5,"mm"))
  draw(typed+typed2+typed3)
 }
