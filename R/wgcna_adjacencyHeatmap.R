#' @title plots adjacency heatmap used by wgcna_analysis
#' @description for analyzing a network we want to plot all the interesting biotypes and look at their individiual correlations and relationship to Module EigeneGenes.
#' @param MEs  module eigen genes
#' @param datTraits data traits
#' @param biotype  all the biotypes of interest
#' @import WGCNA
#' @export
#' @return images
wgcna_adjacencyHeatmap<-function(MEs,datTraits,biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"))
{
 stopifnot(is.null(MEs)==FALSE)

  for(i in 1:length(biotype)){
   weight<-as.data.frame(datTraits[,grep(biotype[i],colnames(datTraits))])
   if(ncol(weight)<1){
   next
   }
   names(weight) = biotype[i]
   MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  plotEigengeneNetworks(MET, paste0("MEs ~ ", biotype[i]), marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
  readkey()

# Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  plotEigengeneNetworks(MET, paste0("Eigengene dendrogram ~",biotype[i]), marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
  readkey()
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(MET, paste0("Eigengene adjacency heatmap ~ ",biotype[i]) , marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
  readkey()
  message(paste0("finished ",biotype[i]))
  
  
}



} #main
