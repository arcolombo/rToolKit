#' @title runs a network analysis at the gene level to investigate repeat initiating pathways 
#' @description running WGNCA could find the gene correlation patterns to infer a relationship behind repeat activation.  what causes repeat activation ?  is it a single mutation?  or protein ?  chromatin accessibility?.  This will create the blockwise/single module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype.
#' @import WGCNA
#' @export
#' @return images and cluster at the gene and repeat level
wgcna<-function(kexp,read.cutoff=2,cutHeight=540000,minBranch=2,whichWGCNA=c("single","block")){
  #FIX ME: cluster a TOM for repeats and coding genes to determine modules.  
  #FIX ME: run enrichment for each module.
  ##prepare data
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  ###rows must be SAMPLES columns genes
  cpm<-collapseBundles(kexp,"gene_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  datExpr0<-t(cpm)
  gsg<-goodSamplesGenes(datExpr0,verbose=3)
  ##rows must be SAMPLES columns repeats
  rexp<-findRepeats(kexp)
  rpm<-collapseBundles(rexp,"tx_biotype",read.cutoff=read.cutoff)
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
  datTraits<-t(rpm)
  datTraits<-as.data.frame(datTraits)
  stopifnot(nrow(datExpr0)==nrow(datTraits))
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

cutHeight<-readHeight()
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
                    main="TxBiotype Correlation Samples") 
    
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

selectedPower<-readPower()
##FIX ME:  add path controls

if(whichWGCNA=="single"){
 ##auto####################################################
net = blockwiseModules(datExpr, power = selectedPower,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "rwaMouseTOM", 
                       verbose = 3)



# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;
}
####compare auto to block
##FIX ME: add hist plot of the many many blocks and colors.
if(whichWGCNA=="block"){
##############BLOCK LEVEL ###################

###first annotate and take out non-annotated first run BWA on datExpr that has annotations, no NAs
message("annotating...")
datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
annot<-geneAnnotation(datExpr,datTraits)
datExpr<-datExpr[,names(datExpr)%in%annot$ensembl_gene_id]
probes2annot<-match(names(datExpr),annot$ensembl_gene_id)
stopifnot(sum(is.na(probes2annot))==0) ##no NA 

##############################################

bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                       power = selectedPower, TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "rwaTOM-blockwise",
                       verbose = 3)
# Load the results of single-block analysis

bwLabels = matchLabels(bwnet$colors,bwnet$colors)
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)



# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# this line corresponds to using an R^2 cut-off of h
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes;
} ##by block

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 11.5, 3, 3));
# Display the correlation values within a heatmap plot

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

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
               main = paste("Module-Repeat Biotype relationships"))
readkey()
colorDF<-sapply(MEs,function(x) median(x))
colorDF<-colorDF[order(colorDF,decreasing=TRUE)]

par(mar = c(11.5, 6, 3, 3));
plot(colorDF,main="Median Correlation Per Module")
axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
readkey()

return(list(datExpr=datExpr,
            datTraits=datTraits,
            medianCor=colorDF,
            annot=annot ))
 
}#main
