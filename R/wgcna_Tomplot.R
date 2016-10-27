#' @title plots topological overlap matrices
#' @description plots the TOM plot given the wgcna matrices
#' @param lnames an object from wgcna method call
#' @param softPower integer 6 or 7 does nice for log2 XR data.
#' @import WGCNA
#' @export
#' @return images of TOM
wgcna_Tomplot<-function(lnames,softPower=6){
 ##FIX ME TOM Plot requires 11.9 gb and will error
 if(is.null(lnames)==TRUE){
  load("wgcna.dataInput.RData")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  geneTree<-lnames[["geneTree"]]
  } else if(is.null(lnames)==FALSE){
  message("found wgcna.dataInput")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  geneTree<-lnames[["geneTree"]]
  } else {
  stop("please run wgcna, need the output before analyzing.")
  }
 
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  message("supporting signed Topological maps and networks only.")
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softPower,
                             corType="bicor",
                             networkType="signed",
                             TOMType="signed" );
  plotTOM= dissTOM^7
  diag(plotTOM)<-NA
  TOMplot(plotTOM, geneTree, bwModuleColors, main = "Network heatmap plot, all genes")  




  adjacency<-adjacency(datExpr,power=softPower)
  TOM<-TOMsimilarity(adjacency)
  dissTOM<-1-TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

} #main
