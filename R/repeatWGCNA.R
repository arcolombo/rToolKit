#' @title runs a network analysis at the repeat level to investigate repeat initiating pathways creating repeat modules agnostically
#' @description running WGNCA could find the repeat correlation patterns to infer a relationship behind repeat activation. This uses the recommended 'biocor' function which is a bi-weight mid-correlation calculation. For normalization, the default uses tmm normalization and a log2 transformation for each cpm for genes, and repeat txBiotypes in the *same* way but on separate calls (this seems ok intuitively).  We use 'signed' networks based on the FAQ.  This will create the blockwise/single module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype.
#' @param kexp a repeat kexp over stages
#' @param read.cutoff integer floor filter
#' @param minBranch integer for cluter min
#' @param whichWGCNA character single or block analysis, block is more sensitive
#' @param entrezOnly boolean, soon to be deprecated because entrez is auto filtered when enrichment testing
#' @param species char, mouse or humans
#' @param selectedPower  6 usually is good. can rerun if NULL
#' @param tmm.norm boolean true
#' @param useBiCor true
#' @import WGCNA
#' @import edgeR
#' @export
#' @return images and cluster at the gene and repeat level
wgcna<-function(kexp,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),entrezOnly=FALSE,species=c("Homo.sapiens","Mus.musculus"),selectedPower=NULL,tmm.norm=TRUE,useBiCor=TRUE,geneMEs=NULL){
  ##FIX ME: add option for TPMs as well as  CPM. Add TPM support
  
  ##prepare data
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  
  ###rows must be SAMPLES columns genes
  cpm<-collapseBundles(kexp,"gene_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
 
 
 if(tmm.norm==TRUE){
  d<-DGEList(counts=cpm)
  cpm.norm<-cpm(d,normalized.lib.sizes=TRUE)
  cpm<-cpm.norm
  cpm.norm<-NULL
   }#tmm norm
  cpm<-log2(1+cpm) ##log2 transform recommended of repeats
  #split out repeats
  datExpr0<-t(cpm)
  gsg<-goodSamplesGenes(datExpr0,verbose=3)
  ##rows must be SAMPLES columns repeats
 
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

 # if((nrow(datExpr)!=nrow(datExpr0))==TRUE){
 #datTraits<-datTraits[keepSamples,]
 # }
# Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")

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
#  plot(sft$fitIndices[,1], sft$fitIndices[,5],
 #    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  #   main = paste("Mean connectivity"))
  selectedPower<-readPower()
  } #selectedPower NULL
  message("annotating...")
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  annot<-as.data.frame(rowRanges(kexp)[rownames(counts(kexp))%in%colnames(datExpr),])
 annot<-data.frame(ensembl_gene_id=annot$tx_id,
                   entrezgene=factor("NA"),
                   hgnc_symbol=annot$tx_biotype,
                   description=annot$gene_biotype)
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
#############################
 enableWGCNAThreads()
net = blockwiseModules(datExpr, power = selectedPower,
                       networkType="signed",
                       corType="bicor",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "rwasingleTOM", 
                       verbose = 3)

# Convert labels to colors for plotting
# Plot the dendrogram and the module colors underneath
  bwLabels<-net$colors ###for saving
  bwModuleColors = labels2colors(net$colors)
  MEs = net$MEs; ##use the module network calculation, do not recalculate 2nd time
  #plots each gene tree one by one
  wgcna_plotAll_dendrograms(bwnet=net,net=net,whichWGCNA="single")
  
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  geneTree = net$dendrograms;
  if(useBiCor==TRUE){
  moduleTraitCor<-bicor(MEs,datTraits)
  moduleTraitPvalue = bicorAndPvalue(MEs,datTraits,use="pairwise.complete.obs",alternative="two.sided")[["p"]]
   modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)
  } else {
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)
  }
 
 lnames <-list(datExpr=datExpr,
            datTraits=datTraits,
            annot=annot,
            MEs=MEs,
            moduleLabels=bwLabels,
            moduleColors=bwModuleColors,
            geneTree=geneTree,
            moduleTraitCor=moduleTraitCor,
            moduleTraitPvalue=moduleTraitPvalue,
            modulePvalFisher=modulePvalFisher,
            usedbiCor=useBiCor)
} ##single block should have 1 module per datTraits column
  if(whichWGCNA=="block"){
##############BLOCK LEVEL ###################
  message("networking...")
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  ##############################################
  
  enableWGCNAThreads()
  bwnet = blockwiseModules(datExpr, 
                       maxBlockSize = 4000,
                       power = selectedPower, 
                       networkType="signed",
                       TOMType = "signed", 
                       corType="bicor",
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
  save(bwnet,file="bwnet.RData",compress=TRUE)
  # open a graphics window
  sizeGrWindow(6,6)
 ########################################################################
  ##plot gene tree one by one 
  wgcna_plotAll_dendrograms(bwnet=bwnet,whichWGCNA="block",bwModuleColors=bwModuleColors)
# this line corresponds to using an R^2 cut-off of h
  # Recalculate MEs with color labels
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  MEs0 = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
  MEs = orderMEs(MEs0)
  if(useBiCor==TRUE){
  moduleTraitCor<-bicor(MEs,datTraits)
  moduleTraitPvalue = bicorAndPvalue(MEs,datTraits,use="all.obs",alternative="two.sided")[["p"]]
  modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)

  } else {
   moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)
   }
  lnames<-list(datExpr=datExpr,
            datTraits=datTraits,
            annot=annot,
            MEs=MEs,
            moduleLabels=bwLabels,
            moduleColors=bwModuleColors,
            geneTree=geneTree,
            moduleTraitCor=moduleTraitCor,
            moduleTraitPvalue=moduleTraitPvalue,
            modulePvalFisher=modulePvalFisher,
            biCor=useBiCor)
  } ##by block

  # Define nmbers of genes and samp
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
  dev.off()
  wgcna_Cormap(lnames,read.cutoff=read.cutoff,plotDot=FALSE) 
  save(lnames,file="wgcna.dataInput.RData",compress=TRUE)
  
  return(lnames)
 
}#main
