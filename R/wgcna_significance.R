#' @title this runs through significance analysis of modules
#' @description this will analyze the significance of modules returning a heatmap, and average gene significance per module
#' @param lnames the output from lnames wgcna call
#' @import WGCNA
#' @export
#' @return significance images
wgcna_significance<-function(lnames,y=NULL,n.modules=7){


colorh1<-lnames[["moduleColors"]]
datExpr<-lnames[["datExpr"]]
datME<-moduleEigengenes(datExpr,colorh1)$eigengenes

  dissimME=(1-t(cor(datME, method="p")))/2
  hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based of the module eigengenes")
  readkey()
    intModules<-findInterestingModules(n.modules=n.modules)

  if(is.null(y)==TRUE){
  y1<-length(grep("pHSC_",rownames(datExpr)))
  y2<-length(grep("LSC_",rownames(datExpr)))
  y3<-length(grep("Blast_",rownames(datExpr)))
  y<-c(rep(1,y1),rep(2,y2),rep(3,y3))
  }
  
  intModules<-paste0("ME",intModules)

sizeGrWindow(8,9)
plotMEpairs(datME[,colnames(datME)%in%intModules ] ,y=y)
##FIX ME: this shows the relationship between modules,  but do this same plot between traits so that you can plot a fixed modules, and compare across repeat types
ModuleEigengeneNetwork<-datME[,colnames(datME)%in%intModules ]
cor.intModules<-signif(cor(datME, ModuleEigengeneNetwork[,-1]),2)
require(ComplexHeatmap)
Heatmap(cor.intModules,
         column_title="Selected Modules Correlation",
          name="correl")
 readkey()

###must call wgcna_heatmap automatically for all intModules %%3
 moduleHeatmap(datExpr,colorh1,intModules)

} #main 
