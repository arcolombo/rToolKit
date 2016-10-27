#' @title plots the stage DE biotype
#' @description plots the stage of each DE
#' @param kexp a kallisot Experiemtn at the stage repeat level
#' @import arkas
#' @import edgeR
#' @export
#' @return a bunch of useful plots
fitTxBiotypes<-function(kexp,design=NULL){
  
   tx_bundles<-collapseBundles(kexp,bundleID="tx_biotype")
  ##fit the bundles
    rps<-list()
    rge<-DGEList(counts=tx_bundles)
    rge<-calcNormFactors(rge)
    rps$design<-design
    rps$voomed<-voom(rge,rps$design)
    rps$fit<-eBayes(lmFit(rps$voomed,rps$design))
   return(rps) 
}
