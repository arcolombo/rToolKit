#' @title produces a correlation map given wgcna network data
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image.
#' @param read.cutoff  integer for min cutoff
#' @import WGCNA
#' @export
#' @return images of eigengenes
wgcna_Cormap<-function(lnames,read.cutoff=2,plotDot=FALSE){

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

# Will display correlations and their p-values
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 11.5, 3, 3));
# Display the correlation values within a heatmap plot
plot.new()
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               yColorWidth=0.07,
               cex.lab.y=.7,
               colors.lab.y=1.3,
               main = paste0("Module-Repeat Biotype TMM cut:",read.cutoff))


readkey()
colorDF<-sapply(MEs,function(x) median(x))
colorDF<-colorDF[order(colorDF,decreasing=TRUE)]
  if(plotDot==TRUE){
  par(mar = c(11.5, 6, 3, 3));
  plot(colorDF,main="Median Correlation Per Module")
  axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
  readkey()
  } else {
  pdf("correlation_plots.pdf")
    par(mar = c(11.5, 6, 3, 3));
    labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
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
