#' @title produces a correlation map given wgcna network data
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image.
#' @param lnames this is the results from wgcna.R
#' @param read.cutoff  integer for min cutoff
#' @param plotDot boolean this plots the median correlation score for each biotypes in a line plot
#' @param recalc boolean if truen then will recalculate the bicor and pvalues
#' @param targets this is a character vector of modules of interest to subset the corMap
#' @param orderBiotype either Alu or L1.
#' @import WGCNA
#' @export
#' @return images of eigengenes
wgcna_Cormap<-function(lnames,read.cutoff=2,plotDot=FALSE,recalc=FALSE,targets=NULL,orderBiotype="Alu"){
  orderBiotype<-match.arg(orderBiotype,c("Alu","L1"))
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
  nGenes= ncol(datExpr);
  nSamples = nrow(datExpr);
 # Recalculate MEs with color labels
  if(recalc==TRUE){
   moduleTraitCor = bicor(MEs, datTraits, use = "all.obs");
   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
   } else{
   moduleTraitCor<-lnames[["moduleTraitCor"]]
   moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
   }
  if(is.null(targets)==FALSE) {
  moduleTraitCor<-moduleTraitCor[rownames(moduleTraitCor)%in%targets,]
  moduleTraitPvalue<-moduleTraitPvalue[rownames(moduleTraitPvalue)%in%targets,]
 }
  if(orderBiotype=="L1"){
 moduleTraitCor<-moduleTraitCor[order(moduleTraitCor[,8],decreasing=TRUE),]
 } else if(orderBiotype=="Alu"){
  moduleTraitCor<-moduleTraitCor[order(moduleTraitCor[,1],decreasing=TRUE),]
 }

 moduleTraitPvalue.id<-match(rownames(moduleTraitCor),rownames(moduleTraitPvalue))
 moduleTraitPvalue<-moduleTraitPvalue[moduleTraitPvalue.id,]
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "");
 
 stopifnot(dim(textMatrix)==dim(moduleTraitCor))
 
 par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
  plot.new()
  labeledHeatmap(Matrix=moduleTraitCor,
,                xLabels=names(datTraits),
                 yLabels= rownames(moduleTraitCor),
                 ySymbols=rownames(moduleTraitCor), 
                 colorLabels=FALSE,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               yColorWidth=0.07,
               cex.lab.y=.7,
               colors.lab.y=1.3,
               main = paste("Module-Repeat Biotype relationships"))


  readkey()
  colorDF<-sapply(MEs,function(x) median(x))
  colorDF<-colorDF[order(colorDF,decreasing=TRUE)]
  if(plotDot==TRUE){
  plot.new()
  par(mar = c(11.5, 6, 3, 3));
  plot(colorDF,main="Median Correlation Per Module")
  axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
  readkey()
  } else {
  pdf("correlation_plots.pdf",width=12,height=12)
    par(mar = c(11.5, 5, 3, 3));
    labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
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
