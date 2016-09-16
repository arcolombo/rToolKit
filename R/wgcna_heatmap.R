#' @title used in wgcna_significance 
#' @description the wgcna_significane has user input to investigate interesting modules. these intModules are saved during this script.  wgcna_heatmap will produce the datExpr values across all samples
#' @param datExpr  a data expression
#' @param colorh1 a color vector labels
#' @param intModules a vector of length 3
#' @import WGCNA
#' @export
#' @return cool heatmap images of modules over samples
moduleHeatmap<-function(datExpr,colorh1,intModules){
 n.calls<-length(intModules)/3
 intModules<-gsub("ME","",intModules)
 for(i in 1:ceiling(n.calls)){
    sizeGrWindow(8,9)
    par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
    start<-3*i-2
    end<-3*i
  for(j in start:end){
    .wgcna_heatmap(datExpr,colorh1,which.module=intModules[j])
       }
  readkey()
  plot.new()
 }
} #main





.wgcna_heatmap<-function(datExpr,colorh1,which.module=NULL){
  plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
         nrgcols=30,
         rlabels=T,
         rcols=which.module,
         title=which.module,
         clabels=rownames(datExpr) )
      

} #wgcna_heat
