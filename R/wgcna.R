#' @title runs a network analysis at the gene level to investigate repeat initiating pathways 
#' @description running WGNCA could find the gene correlation patterns to infer a relationship behind repeat activation. This uses the recommended 'biocor' function which is a bi-weight mid-correlation calculation. For normalization, the default uses tmm normalization and a log2 transformation for each cpm for genes, and repeat txBiotypes as phenotypic data in the *same* way but on separate calls (this seems ok intuitively).  We use 'signed' networks based on the FAQ.  This will create the blockwise/single module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype (phenotypic relationship)
#' @param kexp a kexp 2 group stage is preferred
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
#' @import limma
#' @export
#' @return images and cluster at the gene and repeat level
wgcna<-function(kexp,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),entrezOnly=FALSE,species=c("Homo.sapiens","Mus.musculus"),selectedPower=NULL,intBiotypes=c("acromeric","centromeric","CR1","Alu","DNA transposon","Endogenous Retrovirus","ERV1","ERV3","ERVK","ERVL","hAT","HSFAU","L1","L2","LTR Retrotransposon","Eutr1","Merlin","PiggyBac","Pseudogene","Repetitive element","satellite","snRNA","SVA","TcMar","telo","Transposable Element","Satellite"),useAllBiotypes=FALSE,tmm.norm=TRUE,useBiCor=TRUE,how=c("cpm","tpm"),batchNormalize=FALSE,batchVector=NULL,design=NULL,collapseBy=c("gene_id","gene_name")){
  if(batchNormalize==TRUE && is.null(design)==TRUE){
   stopifnot(is.null(metadata(kexp)$design)==FALSE)
   design<-metadata(kexp)$design
  }
  if(batchNormalize==TRUE && is.null(batchVector)==TRUE){
 stopifnot(is.null(metadata(kexp)$batch)==FALSE)
   batchVector<-metadata(kexp)$batch
 }

  collapseBy=match.arg(collapseBy,c("gene_id","gene_name"))
   how<-match.arg(how,c("cpm","tpm"))
  ##prepare data
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  byWhich<-"gene"  

  if(how=="cpm"){
  ###rows must be SAMPLES columns genes
  if(batchNormalize==FALSE){
  cpm<-collapseBundles(kexp,collapseBy,read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  rexp<-findRepeats(kexp)
  rpm<-collapseBundles(rexp,"tx_biotype",read.cutoff=read.cutoff) 
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
  ##some intBiotypes will be filtered if the rows of rpm don't pass the read.cutoff
 
 if(tmm.norm==TRUE){
  d<-DGEList(counts=cpm)
  d<-calcNormFactors(d)
  cpm.norm<-cpm(d,normalized.lib.sizes=TRUE,log=FALSE)
  cpm<-cpm.norm
  cpm.norm<-NULL
  rd<-DGEList(counts=rpm)
  rd<-calcNormFactors(rd)
  rdm.norm<-cpm(rd,normalized.lib.sizes=TRUE,log=FALSE)
  rpm<-rdm.norm
  rdm.norm<-NULL
      }#tmm norm
    rpm<-log2(1+rpm)
    intBiotypes<-rownames(rpm)
    cpm<-log2(1+cpm)
    }else if(batchNormalize==TRUE){
   ##task: input kexp, metadata$design, metadata$batch
   ##the caller will batch normalize WITH design matrix
   ##output: log2 batch cpm, re-convert into CPM, collapse tx_id normalized into gene bundles CPM ,final output is batch correct CPM values tx_id
  ##taks for this branch is to collapse tx_id into gene_id (CPM), and tx_biotype(cpm)
   cpm<-collapseBundles(kexp,"tx_id",read.cutoff=read.cutoff)
   d<-DGEList(counts=cpm)
   d<-calcNormFactors(d)
   log.cpm<-cpm(d,log=TRUE,normalize.lib.sizes=TRUE) #log2
    ##removeBatch requires log2
   logCPM<-removeBatchEffect(log.cpm,batch=batchVector,design=design)
   CPM<-logCPM^2
   ####now split gene and repeats
 
    bundleable <- !is.na(mcols(rowRanges(kexp))[["gene_id"]])
    feats<-rowRanges(kexp)[bundleable] 
    feats<-feats[rownames(CPM),]
    bundleable<-!is.na(mcols(feats)[["gene_id"]])
  cts <- split.data.frame(CPM[bundleable, ], mcols(feats)[["gene_id"]])
  cts <- cts[ sapply(cts, function(x) max(x) >= read.cutoff) ]
  cts <- lapply(cts, function(x) x[ rowSums(x) >= read.cutoff, ])
 bundled <- do.call(rbind, lapply(cts,
                                   function(x)
                                     if(!is.null(nrow(x))) colSums(x) else x))
   ####split genes and repeats
  cpm.norm<-bundled[grepl("^ENS",rownames(bundled)),]
  rpm.norm<-bundled[!grepl("^ENS",rownames(bundled)),]
  rpm.norm<-rpm.norm[!grepl("^ERCC",rownames(rpm.norm)),] ##repeat batch free tx
  ##########collapse batch free repeat transcripts into tx_biotype aggregates
    bundleable <- !is.na(mcols(rowRanges(kexp))[["tx_biotype"]])
    feats<-rowRanges(kexp)[bundleable]
    feats<-feats[rownames(rpm.norm),]
    bundleable<-!is.na(mcols(feats)[["tx_biotype"]])
  rts <- split.data.frame(rpm.norm[bundleable, ], mcols(feats)[["tx_biotype"]])
  rts <- rts[ sapply(rts, function(x) max(x) >= read.cutoff) ]
  rts <- lapply(rts, function(x) x[ rowSums(x) >= read.cutoff, ])
 tx.bundled <- do.call(rbind, lapply(rts,
                                   function(x)
                                     if(!is.null(nrow(x))) colSums(x) else x))

  cpm<-log2(1+cpm.norm)
  rpm<-log2(1+tx.bundled)
  ############collapse by gene_id
  }##batch
  }else if(how=="tpm"){
  cpm<-collapseTpm(kexp,collapseBy,read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  rexp<-findRepeats(kexp)
  rpm<-collapseTpm(rexp,"tx_biotype",read.cutoff=read.cutoff)
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
  cpm<-log2(1+cpm) ##log2 transform recommended of genes
  rpm<-log2(1+rpm) ##log2 transform of repeats.
  }
  intBiotypes<-intBiotypes[intBiotypes%in%rownames(rpm)]
  #split out repeats
  datExpr0<-t(cpm)
  gsg<-goodSamplesGenes(datExpr0,verbose=3)
  ##rows must be SAMPLES columns repeats
 
  if(useAllBiotypes==FALSE){
  #select columns of intBiotypes
   stopifnot(is.null(intBiotypes)==FALSE)
   rpm<-rpm[rownames(rpm)%in%intBiotypes,]
    }
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

  pdf(paste0("hclust_",selectedPower,"_",how,".pdf"),width=12,height=9)
  plot(sampleTree, main = "Sample clustering to detect outliers",
       sub="",
      xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
   dev.off()
 

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
   if(useAllBiotypes==FALSE){
   stopifnot(all(intBiotypes%in%colnames(datTraits))==TRUE)
   datTraits<-datTraits[,match(intBiotypes,colnames(datTraits))]
   traitColors<-numbers2colors(datTraits[,match(intBiotypes,colnames(datTraits))],signed=FALSE)
   }else{
   traitColors = numbers2colors(datTraits, signed = FALSE);
   }
# Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    marAll=c(1,11,3,3),
                    main=paste0("TxBiotype ",how," Correlation Samples")) 
  readkey()   
  pdf(paste0("TxBiotype_",how,"_Correlation_Samples.pdf"),width=12,heigh=10)
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    marAll=c(1,11,3,3),
                    main=paste0("TxBiotype ",how," Correlation Samples") )
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
   text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers,cex=cex1,col="red");
  selectedPower<-readPower()
  ####Print to PDF############
  pdf(paste0("WGCNAselectedPower_",how,"analysis.pdf"))
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
   plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
   text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers,cex=cex1,col="red");
  dev.off()

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
   wgcna_plotAll_dendrograms(bwnet=net,whichWGCNA="single",bwModuleColors=bwModuleColors,bwLabels=bwLabels,how=how,byWhich=byWhich)
  
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
            usedbiCor=useBiCor,
            how=how,
            byWhich="gene")
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
  save(bwnet,file=paste0("bwnet_",selectedPower,".RData"),compress=TRUE)
  # open a graphics window
  sizeGrWindow(6,6)
 ########################################################################
  ##plot gene tree one by one 
  wgcna_plotAll_dendrograms(bwnet=bwnet,whichWGCNA="block",bwModuleColors=bwModuleColors,bwLabels=bwLabels,how=how,byWhich="gene")
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
            biCor=useBiCor,
            how=how, 
            byWhich="gene")
  } ##by block
  save(lnames,file=paste0("wgcna.",how,"_",selectedPower,".dataInput.RData"),compress=TRUE)
  dev.off()  
# Display the correlation values within a heatmap plot
  #wgcna_Cormap(lnames,read.cutoff=read.cutoff,plotDot=FALSE,how=how) 
  return(lnames)
 
}#main
