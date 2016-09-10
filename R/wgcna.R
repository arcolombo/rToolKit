#' @title runs a network analysis at the gene level to investigate repeat initiating pathways 
#' @description running WGNCA could find the gene correlation patterns to infer a relationship behind repeat activation.  what causes repeat activation ?  is it a single mutation?  or protein ?  chromatin accessibility?
#' @import WGCNA
#' @export
#' @return images and cluster at the gene and repeat level
wgcna<-function(kexp,read.cutoff=1,cutHeight=540000,minBranch=2){
  #FIX ME: cluster a TOM for repeats and coding genes to determine modules.  
  #FIX ME: run enrichment for each module.
  ##prepare data

  ###rows must be SAMPLES columns genes
  cpm<-collapseBundles(kexp,"gene_id",read.cutoff=read.cutoff)
  datExpr0<-t(cpm)
  gsg<-goodSamplesGenes(datExpr0,verbose=3)
  ##rows must be SAMPLES columns repeats
  rexp<-findRepeats(kexp)
  rpm<-collapseBundles(rexp,"tx_id",read.cutoff=read.cutoff)
  datTraits<-t(rpm)
  datTraits<-as.data.frame(datTraits)
  stopifnot(nrow(datExpr0)==nrow(traitData))
  ##the rows of each clinical/trait/repeat data must match both objects
  
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = cutHeight, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = minBranch)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

if((nrow(datExpr)!=nrow(datExpr0))==TRUE){
 datTraits<-datTraits[keepSamples,]

}
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                     marAll=c(1,11,3,3),
                     main = "Sample dendrogram and trait heatmap")




}#main
