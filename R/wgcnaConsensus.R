#' @title runs a consensus network analysis at the gene level to investigate repeat initiating pathways 
#' @description runs WGNCA consensus analysis could find the gene correlation patterns to infer a relationship behind repeat activation. This uses the recommended 'biocor' function which is a bi-weight mid-correlation calculation. For normalization, the default uses tmm normalization and a log2 transformation for each cpm for genes, and repeat txBiotypes in the *same* way but on separate calls (this seems ok intuitively).  We use 'signed' networks based on the FAQ.  This will create the blockwise/single consensus module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype. Note this inputs 2 kallistoExperiments and creates a multiExpression structure as a necessary pre-processing step for blockWiseConsensus network call. 
#' @param kexp1 a kexp one group stage is preferred
#' @param kexp2 a kexp for consensus analysis
#' @param read.cutoff integer floor filter
#' @param minBranch integer for cluter min
#' @param whichWGCNA character single or block analysis, block is more sensitive
#' @param entrezOnly boolean, soon to be deprecated because entrez is auto filtered when enrichment testing
#' @param species char, mouse or humans
#' @param selectedPower  6 usually is good. can rerun if NULL
#' @param intBiotypes character the tx_biotypes of interest
#' @param useAllBiotypes boolean if false then intBiotypes are used, if true than the correlations are checked against all tx_biotypes
#' @import WGCNA
#' @import edgeR
#' @export
#' @return images and cluster at the gene and repeat level
wgcnaConsensus<-function(kexp1,kexp2,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),entrezOnly=FALSE,species=c("Homo.sapiens","Mus.musculus"),selectedPower=NULL,intBiotypes=c("Alu","DNA transposon","Endogenous Retrovirus","ERV1","ERV3","ERVK","ERVL","L1","L2","LTR Retrotransposon","Satellite"),useAllBiotypes=FALSE,tmm.norm=TRUE,useBiCor=TRUE,setLabels=c("kexp1","kexp2")){
  ##FIX ME: add option for TPMs as well as  CPM. Add TPM support
  
  ##prepare kexp1 and kexp2 for multiExpr data list
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  
  shortLabels<-c("Buenrostro","TCGA")
  ###rows must be SAMPLES columns genes
  cpm1<-collapseBundles(kexp1,"gene_id",read.cutoff=read.cutoff)
  cpm1<-cpm1[!grepl("^ERCC",rownames(cpm1)),]
  cpm2<-collapseBundles(kexp2,"gene_id",read.cutoff=read.cutoff)
  cpm2<-cpm2[!grepl("^ERCC",rownames(cpm2)),]
  ##for a consensus network analysis to work you can only examine common genes between kexp 1 and kexp2
  if(nrow(cpm1)<=nrow(cpm2)){
  id<-match(rownames(cpm1),rownames(cpm2))
  cpm2<-cpm2[id,]
  } else{
    id<-match(rownames(cpm1),rownames(cpm2))
   cpm1<-cpm1[id,] 
 }
  stopifnot(all(rownames(cpm1)%in%rownames(cpm2))==TRUE)
####FIX ME: not sure how to deal with two different trait data sets.
 # rexp1<-findRepeats(kexp1)
 # rpm<-collapseBundles(rexp,"tx_biotype",read.cutoff=read.cutoff) 
 # rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
 if(tmm.norm==TRUE){
  d1<-DGEList(counts=cpm1)
  cpm.norm1<-cpm(d1,normalized.lib.sizes=TRUE)
  cpm1<-cpm.norm1
  cpm.norm1<-NULL
  #rd<-DGEList(counts=rpm)
  #rdm.norm<-cpm(rd,normalized.lib.sizes=TRUE)
  #rpm<-rdm.norm
  #rdm.norm<-NULL
  d2<-DGEList(counts=cpm2)
  cpm.norm2<-cpm(d2,normalized.lib.sizes=TRUE)
  cpm2<-cpm.norm2
  cpm.norm2<-NULL
  }#tmm norm
  cpm1<-log2(1+cpm1) ##log2 transform recommended of genes kexp1
  cpm2<-log2(1+cpm2) ##log2 transform of genes kexp2
  
  # Form multi-set expression data: columns starting from 9 contain actual expression data.
  multiExpr = vector(mode = "list", length = 2)
  multiExpr[[1]] = list(data = as.data.frame(t(cpm1)));
  names(multiExpr[[1]]$data) = rownames(cpm1);
  rownames(multiExpr[[1]]$data) = colnames(cpm1)   ;
  multiExpr[[2]] = list(data = as.data.frame(t(cpm2)));
  names(multiExpr[[2]]$data) = rownames(cpm2);
  rownames(multiExpr[[2]]$data) = colnames(cpm2);
  # Check that the data has the correct format for many functions operating on multiple sets:
  exprSize = checkSets(multiExpr)
  print(exprSize)
  
  gsg<-goodSamplesGenesMS(multiExpr,verbose=3)
  ##rows must be SAMPLES columns repeats
  

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}
###########

  if(useAllBiotypes==FALSE){
  #select columns of intBiotypes
   stopifnot(is.null(intBiotypes)==FALSE)
   rpm<-rpm[rownames(rpm)%in%intBiotypes,]
   
   }
  datTraits<-t(rpm)
  datTraits<-as.data.frame(datTraits)
  stopifnot(nrow(datExpr0)==nrow(datTraits))
  ##the rows of each clinical/trait/repeat data must match both objects
  




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
  pdf("TxBiotype_Correlation_Samples.pdf")
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    marAll=c(1,11,3,3),
                    main="TxBiotype Correlation Samples") 
  dev.off() 
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
  wgcna_plotAll_dendrograms(net,whichWGCNA="single")
  
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
