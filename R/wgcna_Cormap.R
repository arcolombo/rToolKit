#' @title produces a correlation map given wgcna network data
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image.
#' @param read.cutoff  integer for min cutoff
#' @import WGCNA
#' @export
#' @return images of eigengenes
wgcna_Cormap<-function(lnames,read.cutoff=2,plotDot=FALSE,recalc=FALSE){

if(is.null(lnames)==TRUE){
load("wgcna.dataInput.RData")
}
###declare needed objects from load
if(is.null(lnames)==FALSE){
message(paste0("found lnames"))
} else {
 stop("did not find lnames, please run wgnca")
}
bwModuleColors<-lnames[["moduleColors"]]
MEs<-lnames[["MEs"]]
datExpr<-lnames[["datExpr"]]
datTraits<-lnames[["datTraits"]]
annot<-lnames[["annot"]]
moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
# Will display correlations and their p-values
# Define numbers of genes and samples
nGenes= ncol(datExpr);
tmes<-rownames(moduleTraitCor[order(moduleTraitCor[,1],decreasing=TRUE),])
extMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               yColorWidth=0.07,
               cex.lab.y=.7,
               colors.lab.y=1.3,
               main = paste("Module-Repeat Biotype relationships"))
  plot(colorDF,main="Median Correlation Per Module")
  axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
  dev.off()
 }
}
