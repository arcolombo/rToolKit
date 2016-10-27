#' @title visualizes barplot of each module across samples
#' @description The first PC is referred to as the module eigengene (ME), and is a single value that represents the highest percent of variance for all genes in a module. In other words, if we show the ME for module X doing something, there is a good chance that most genes in module X also do that same thing. 
#' @export
moduleVisualization<-function(pcTree=NULL,MDS_G=NULL,colorsG=NULL,treeOrder=NULL,modulesG=NULL,moduleEigen=NULL,imageName="G1"   ){

  par(mfrow=c(1,2), mar=c(0, 3, 1, 1) + 0.1, cex=1)
  plot(pcTree, xlab="",ylab="",main="",sub="")
  plot(MDS_G, col= colorsG, main="MDS plot", cex=2, pch=19)
  readkey()
   ordergenes = treeOrder
   ##skip heatmap, show module expression as bar plot.
   for (which.module in names(table(modulesG))){
   ME = moduleEigen[, paste("ME",which.module, sep="")]
   barplot(ME,
          col=which.module,
          cex.main=2,
          ylab="eigengene expression",
          xlab="array sample",
           main=paste0(names(table(modulesG)[which.module])))
       readkey()
       };


  pdf(paste0(imageName,"_","ModuleEigengeneVisualizations.pdf")
            ,height=6,width=6)
  par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
  plot(pcTree, xlab="",ylab="",main="",sub="")
  plot(MDS_G, col= colorsG, main="MDS plot", cex=2, pch=19)
  ordergenes =treeOrder
   ##skip heatmap, show module expression as bar plot.
   for (which.module in names(table(modulesG))){
   ME = moduleEigen[, paste("ME",which.module, sep="")]
   barplot(ME,
          col=which.module,
          cex.main=2,
          ylab="eigengene expression",
          xlab="array sample",
           main=paste0(names(table(modulesG)[which.module])))
           };
   dev.off()
 }
