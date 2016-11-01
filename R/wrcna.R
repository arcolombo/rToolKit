#' @title creates a network module at the repeat level to investigate repeat initiating pathways 
#' @description this differs from wgcna.R, where wgcna.R create phenotypic module analysis of genes and phenotypic repeat data.  wrna will run a single/block wise adjacency/correlation matrices of repeats into modules of transcripts and also returns an object ready for wgcnaDbLite that includes repeat phenotypic data. This uses the recommended 'biocor' function which is a bi-weight mid-correlation calculation. For normalization, the default uses tmm normalization and a log2 transformation for each cpm for genes,  We use 'signed' networks based on the FAQ.  This will create the blockwise/single module data frame and run the soft-thresholding and create a correlation heatmap.  the downstream method is wgcna_analsyis which investigates specific module color and specific biotype (phenotypic relationship).  tx_biotype is included as phenotypic data
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
#' @export
#' @return images and cluster at the gene and repeat level
wrcna<-function(kexp,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),species=c("Homo.sapiens","Mus.musculus"),selectedPower=6,intBiotypes=c("Alu","DNA transposon","Endogenous Retrovirus","ERV1","ERV3","ERVK","ERVL","L1","L2","LTR Retrotransposon","Satellite"),useAllBiotypes=FALSE,tmm.norm=TRUE,useBiCor=TRUE,how=c("cpm","tpm"), batchNormalize=FALSE,batchVector=NULL){
  
   how<-match.arg(how,c("cpm","tpm"))
   byWhich<-"repeat"
  ##prepare data
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  rexp<-findRepeats(kexp)


  if(how=="cpm"){
  cpm<-collapseBundles(rexp,"tx_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  cpm<-cpm[!grepl("^ENS",rownames(cpm)),]
  rpm<-collapseBundles(rexp,"tx_biotype",read.cutoff=read.cutoff) 
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
   if(tmm.norm==TRUE){
  d<-DGEList(counts=cpm)
  cpm.norm<-cpm(d,normalized.lib.sizes=TRUE)
  cpm<-cpm.norm
  cpm.norm<-NULL
  rd<-DGEList(counts=rpm)
  rdm.norm<-cpm(rd,normalized.lib.sizes=TRUE)
  rpm<-rdm.norm
  rdm.norm<-NULL
   if(batchNormalize==TRUE){
     stopifnot(is.null(batchVector)==FALSE)
     stopifnot(length(batchVector)==ncol(kexp))##the batchVector nomenclature must match the length of columns.
     ##takes TMM normalized and batch corrects
     batch.cpm<-removeBatchEffect(cpm,batch=batchVector)
     batch.rpm<-removeBatchEffect(rpm,batch=batchVector)
      }
     }#tmm.norm
  }else if(how=="tpm"){
  cpm<-collapseTpm(rexp,"gene_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  cpm<-cpm[!grepl("^ENS",rownames(cpm)),]
  rpm<-collapseTpm(rexp,"tx_biotype",read.cutoff=read.cutoff)
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
  }
  if(batchNormalize==FALSE){
  cpm<-log2(1+cpm) 
  rpm<-log2(1+rpm) ##log2 transform of repeats.
  }else if(batchNormalize==TRUE){
   cpm<-log2(1+batch.cpm)
    rpm<-log2(1+batch.rpm)
    ### filter NaNs
    id<-which(is.na(cpm))
    cpm[id]<-1.0
    rid<-which(is.na(rpm))
    rpm[rid]<-1.0   
    }
  datExpr0<-t(cpm)

   gsg<-goodSamplesGenes(datExpr0,verbose=3)
  ##rows must be SAMPLES columns repeats
  
  if(useAllBiotypes==FALSE&&byWhich=="repeat"){
  #select columns of intBiotypes
   stopifnot(is.null(intBiotypes)==FALSE)
   rpm<-rpm[rownames(rpm)%in%intBiotypes,]
    }
  datTraits<-t(rpm)
  datTraits<-as.data.frame(datTraits)
  stopifnot(nrow(datExpr0)==nrow(datTraits))

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
  plot(sampleTree, main =paste0("Sample ",byWhich," clustering to detect outliers"), 
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
  print(paste0("samples to omit ",colnames(rexp)[which(keepSamples==FALSE)]))
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
   if((nrow(datExpr)!=nrow(datExpr0))==TRUE){
 datTraits<-datTraits[keepSamples,]
  }
 
# Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
   traitColors = numbers2colors(datTraits, signed = FALSE);

  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    marAll=c(1,11,3,3),
                    main=paste0("Repeat ",how," Module TxBiotype Correlation Samples"))
  readkey()
  pdf(paste0("RepeatMM_",how,"_TxBiotype_Correlation_Samples.pdf"),width=12,height=9)
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    marAll=c(1,11,3,3),
                    main=paste0("Repeat ",how," Module TxBiotype Correlation Samples"))
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
  y<-( -sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  x<-sft$fitIndices[,1] 
 plot(x,y,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type='n',  
   main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  readkey()
   y2<-sft$fitIndices[,5]
  plot(x, y2,
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
   text(x, y2,
     labels=powers,cex=cex1,col="red");

  selectedPower<-readPower()
  pdf(paste0("RepeatModule_",how,"_soft_ThresholdPower.pdf"),width=12,height=9)
  plot(x,y,
     xlab=paste0("RE ",how," Soft Threshold (power)"),
     ylab="RE Scale Free Topology Model Fit,signed R^2",
     type='n',
   main = paste("RE Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
    abline(h=0.80,col="red")
    y2<-sft$fitIndices[,5]
  plot(x, y2,
     xlab="RE Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("RE Mean connectivity"))
   text(x, y2,
     labels=powers,cex=cex1,col="red");
   dev.off()

  } #selectedPower NULL
  message("annotating...")
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  ####### ensembl_gene_id  entrezgene hgnc_symbol description   add data here
  annot<-DataFrame(rowRanges(rexp)[!duplicated(rowRanges(rexp)$tx_id)])
  txID<-grep("tx_id",colnames(annot))
  entrezID<-grep("entrezid",colnames(annot))
  geneID<-grep("gene_id",colnames(annot))
  txBioID<-grep("tx_biotype",colnames(annot))
  annot<-annot[,c(txID,entrezID,geneID,txBioID)]
  annot<-annot[,!grepl("^X.",colnames(annot))]
  colnames(annot)<-c("ensembl_gene_id","entrezid","hgnc_symbol","description")
  annot$entrezgene<-"NA"
  annot<-as.data.frame(annot)



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
  
 rnames <-list(datExpr=datExpr,
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
            byWhich=byWhich)
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
  wgcna_plotAll_dendrograms(bwnet=bwnet,whichWGCNA="block",bwModuleColors=bwModuleColors,bwLabels=bwLabels,how=how,byWhich=byWhich)
# this line corresponds to using an R^2 cut-off of h
  # Recalculate MEs with color labels
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  MEs0 = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
  MEs = orderMEs(MEs0)
  if(useBiCor==TRUE){
   moduleTraitCor<-bicor(MEs,datTraits)
  moduleTraitPvalue = bicorAndPvalue(MEs,datTraits,use="pairwise.complete.obs",alternative="two.sided")[["p"]]
    modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)
  } else {
   moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  modulePvalFisher<-corPvalueFisher(moduleTraitCor,nSamples)
  }


  rnames<-list(datExpr=datExpr,
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
            byWhich=byWhich)
  } ##by block
  save(rnames,file=paste0("wgcna.",how,"_",byWhich,".dataInput.RData"),compress=TRUE)
   cat("done.\n")
   dev.off()
   return(rnames)
 }#main
