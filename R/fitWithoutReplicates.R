#' @title this analysis dendritic cells from Zhang et. al
#' @description dendritic tet2 cells in a 2 group study. this analyzes 2 groups 1 comparison and 1 control. this performs an MDS plot of the dendritic 4h/0h 
#' @param kexp a kexp of dendritic 2 group 
#' @param bundleID  collapse range column
#' @param read.cutoff  floor expression
#' @export
#' @return a DGE list
#' @import limma
#' @import edgeR
fitWithoutReplicates<-function(kexp,bundleID="tx_id",read.cutoff=1,byWhat=c("genes","repeats"),how=c("counts","tpm"),...){
##FIX ME support low input
  res<-list()
  kexp<-fixRepeatRanges(kexp)
  if(byWhat=="repeat"){
  kexp<-findRepeats(kexp)
  }
  if(how=="counts"){
  bundledCounts<-collapseBundles(kexp,bundleID=bundleID,
                                read.cutoff=read.cutoff,...)
  } else {
  bundledCounts<-collapseTpm(kexp,bundleID="tx_id")
  }
 
  dge<-DGEList(counts=bundledCounts)
  dge<-calcNormFactors(dge)
  plotMA(dge)
  return(dge)
}
