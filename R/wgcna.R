#' @title runs a network analysis at the gene level to investigate repeat initiating pathways 
#' @description running WGNCA could find the gene correlation patterns to infer a relationship behind repeat activation.  what causes repeat activation ?  is it a single mutation?  or protein ?  chromatin accessibility?.  This will create the blockwise/single module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype.
#' @param kexp a kexp 2 group stage is preferred
#' @param read.cutoff integer floor filter
#' @param minBranch integer for cluter min
#' @param whichWGCNA character single or block analysis, block is more sensitive
#' @param entrezOnly boolean, soon to be deprecated because entrez is auto filtered when enrichment testing
#' @param species char, mouse or humans
#' @param selectedPower  6 usually is good. can rerun if NULL
#' @import WGCNA
#' @export
#' @return images and cluster at the gene and repeat level
wgcna<-function(kexp,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),entrezOnly=FALSE,species=c("Homo.sapiens","Mus.musculus"),selectedPower=6){
  ##FIX ME: use TPMs instead of CPM???
  
  ##prepare data
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))

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
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub="", 
      xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

  cutHeight<-readHeight()
# Plot a line to show the cut
  abline(h = cutHeight, col = "red");
  readkey()
# Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = minBranch)
  print(table(clust))

# clust 1 contains the samples we want to keep.
  keepSamples = (clust!=0)
  print(paste0("samples to omit ",colnames(kexp)[which(keepSamples==FALSE)]))
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
readkey()    
# Choose a set of soft-thresholding powers
  if(is.null(selectedPower)==TRUE){
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
  readkey()
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
  selectedPower<-readPower()
  } #selectedPower NULL
  message("annotating...")
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  annot<-geneAnnotation(datExpr,datTraits,species=species)

  if(entrezOnly==TRUE){
  id<-!is.na(annot$entrezgene)
  datExpr<-datExpr[,names(datExpr)%in%annot$ensembl_gene_id[id]]
  probes2annot<-match(names(datExpr),annot$ensembl_gene_id)
  stopifnot(sum(is.na(probes2annot))==0) ##no NA 
  } else{
  datExpr<-datExpr[,names(datExpr)%in%annot$ensembl_gene_id]
  probes2annot<-match(names(datExpr),annot$ensembl_gene_id)
  stopifnot(sum(is.na(probes2annot))==0) ##no NA 
  }


##FIX ME:  add path controls

if(whichWGCNA=="single"){
 ##auto####################################################
datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)



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
bwLabes<-net$colors ###for saving
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
bwModuleColors<-moduleColors ##for saving
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;

####compare auto to block
##FIX M
lnames<-list(datExpr=datExpr,
            datTraits=datTraits,
            medianCor=colorDF,
            annot=annot,
            MEs=MEs,
            moduleLabels=bwLabels,
            moduleColors=bwModuleColors,
            geneTree=geneTree)
}
  if(whichWGCNA=="block"){
##############BLOCK LEVEL ###################
  message("networking...")
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  ##############################################

  bwnet = blockwiseModules(datExpr, 
                       maxBlockSize = 2000,
                       power = selectedPower, 
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "rwaTOM-blockwise",
                       verbose = 3)
# Load the results of single-block analysis

  bwLabels = matchLabels(bwnet$colors,bwnet$colors)
  # Convert labels to colors for plotting
  bwModuleColors = labels2colors(bwLabels)
  geneTree<-bwnet$dendrograms
  moduleColors<-bwModuleColors

  # open a graphics window
  sizeGrWindow(6,6)
  # Plot the dendrogram and the module colors underneath for block 1
  plotDendroAndColors(bwnet$dendrograms[[1]], 
                      bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", 
                     main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
  plotDendroAndColors(bwnet$dendrograms[[2]], 
                      bwModuleColors[bwnet$blockGenes[[2]]],
                      "Module colors", 
                      main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# this line corresponds to using an R^2 cut-off of h
  blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes;
  } ##by block

  sizeGrWindow(10,6)
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
  lnames<-list(datExpr=datExpr,
            datTraits=datTraits,
            medianCor=colorDF,
            annot=annot,
            MEs=MEs,
            moduleLabels=bwLabels,
            moduleColors=bwModuleColors,
            geneTree=geneTree )
  save(lnames,file="wgcna.dataInput.RData",compress=TRUE)
  return(list(datExpr=datExpr,
            datTraits=datTraits,
            medianCor=colorDF,
            annot=annot,
            MEs=MEs,
            moduleLabels=bwLabels,
            moduleColors=bwModuleColors,
            geneTree=geneTree ))
 
}#main
