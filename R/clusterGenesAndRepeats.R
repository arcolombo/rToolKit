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
clusterGenesAndRepeats<-function(kexp,read.cutoff=2,minBranch=2,whichWGCNA=c("single","block"),entrezOnly=FALSE,species=c("Homo.sapiens","Mus.musculus"),selectedPower=NULL,intBiotypes=c("acromeric","centromeric","CR1","Alu","DNA transposon","Endogenous Retrovirus","ERV1","ERV3","ERVK","ERVL","hAT","HSFAU","L1","L2","LTR Retrotransposon","Eutr1","Merlin","PiggyBac","Pseudogene","Repetitive element","satellite","snRNA","SVA","TcMar","telo","Transposable Element","Satellite"),useAllBiotypes=FALSE,tmm.norm=TRUE,useBiCor=TRUE,how=c("cpm","tpm"),batchNormalize=FALSE,batchVector=NULL,design=NULL,collapseBy=c("gene_id","gene_name")){
 
 ##FIX ME: use numbersToColors on a sorted biotype column to get a single bar plot for a legened and print that out separately.  ##
 if(batchNormalize==TRUE && is.null(design)==TRUE){
   stopifnot(is.null(metadata(kexp)$design)==FALSE)
   design<-metadata(kexp)$design
  }
  if(batchNormalize==TRUE && is.null(batchVector)==TRUE){
 stopifnot(is.null(metadata(kexp)$batch)==FALSE)
   batchVector<-metadata(kexp)$batch
 }
   how<-match.arg(how,c("cpm","tpm"))
   collapseBy<-match.arg(collapseBy,c("tx_id","gene_id","gene_name"))
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
    cpm<-log2(1+cpm)
    }else if(batchNormalize==TRUE){
   ##task: input kexp, metadata$design, metadata$batch
   ##the caller will batch normalize WITH design matrix
   ##output: log2 batch cpm, re-convert into CPM, collapse tx_id normalized into gene bundles CPM ,final output is batch correct CPM values tx_id
  ##taks for this branch is to collapse tx_id into gene_id (CPM), and tx_biotype(cpm)
   cpm<-collapseBundles(kexp,collapseBy,read.cutoff=read.cutoff)
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


datExpr<-datExpr0
# Determine cluster under the line
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
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
                    main=paste0("Hierarchial Cluster of Gene Expression (CPM)"),ylab="Average Euclidean Distance Gene Expression")
  readkey()
  pdf(paste0("TxBiotype_",how,"_Correlation_Samples_",collapseBy,".pdf"),width=12,heigh=10)
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    marAll=c(1,11,3,3),ylab=NULL,main=NULL)
  dev.off()
}
