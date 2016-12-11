#' @title uses limma to normalize data
#' @description this is used by qusageTablesFromWGCNA to call batch normalized enrichment.  this returns two levels of batch free log2 cpm gene_id or gene_name. for qusageTablesFromWGCNA calls the correct level is gene_name.  gene_id can be used in wgcna/wrcna however not currently implemented.  if collapseLevel is gene_id then this outputs a list of batch free gene_id log2 CPM and batch free tx_biotype log2 CPM
#' @param kexp a kallisto experiment 
#' @param batchVector batch factor
#' @param design design matrix of treatment groups without batch covariate
#' @param read.cutoff floor threshold
#' 
#' @param byLevel gene currently supported
#' @import edgeR
#' @import arkas
#' @import limma
#' @return list of batch free log2 cpm and batch free repeat biotypes log2 cpm
batchEffectNormalize<-function(kexp,batchVector=NULL,design=NULL,read.cutoff=2,collapseByLevel=c("gene_id","gene_name")){
 ##generalize this 
###task: intput and batch correct tx_id matrix counts
### input: TMM log2 CPM -> removeBatchEffect
##task :  re-convert batch free log2 ---> CPM batch free (tx_id)
####then collapse batch free CPM into gene_id, and tx_biotype (Perhaps another method, or maybe arkas can use the kexp information to properly generalize a collapse call on a matrix using a kexp.   we can collapse a kexp, but how to collapse a matrix? 

collapseByLevel<-match.arg(collapseByLevel,c("gene_id","gene_name"))
 if(batchNormalize==TRUE && is.null(design)==TRUE){
   stopifnot(is.null(metadata(kexp)$design)==FALSE)
   design<-metadata(kexp)$design
  }
  if(batchNormalize==TRUE && is.null(batchVector)==TRUE){
 stopifnot(is.null(metadata(kexp)$batch)==FALSE)
   batchVector<-metadata(kexp)$batch
 }

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
     ###CPM is in tx_id form which matches feats
    feats<-feats[rownames(CPM),]
    bundleable<-!is.na(mcols(feats)[["gene_id"]])
  cts <- split.data.frame(CPM[bundleable, ], mcols(feats)[[collapseByLevel]])
  cts <- cts[ sapply(cts, function(x) max(x) >= read.cutoff) ]
  cts <- lapply(cts, function(x) x[ rowSums(x) >= read.cutoff, ])
 bundled <- do.call(rbind, lapply(cts,
                                   function(x)
                                     if(!is.null(nrow(x))) colSums(x) else x))
 if(collapseByLevel=="gene_id"){
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
  batch.list<-list(cpm=cpm,rpm=rpm)
 }else if(collapseByLevel=="gene_name"){
   cpm.norm<-log2(1+bundled)
  batch.list<-list(cpm=cpm.norm) 
  } 

 return(batch.list )

}
