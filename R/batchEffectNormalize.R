#' @title uses combat to normalize data
#' @description can be called to use edgeR batchEffect removal
#' @import edgeR
#' @import arkas
#' @return batch cpm log2 transformed. and a pdf to show results.
batchEffectNormalize<-function(kexp,batch=NULL,design=NULL,read.cutoff=2,byWhich=c("gene_id","tx_id"),byLevel=c("gene","repeat")){
 ##generalize this 
####task: intput and batch correct tx_id matrix counts
#### input: TMM log2 CPM -> removeBatchEffect
###task :  re-convert batch free log2 ---> CPM batch free (tx_id)
#####then collapse batch free CPM into gene_id, and tx_biotype (Perhaps another method, or maybe arkas can use the kexp information to properly generalize a collapse call on a matrix using a kexp.   we can collapse a kexp, but how to collapse a matrix? 

  byLevel<-match.arg(byLevel,c("gene","repeat"))
  byWhich<-match.arg(byWhich,c("gene_id","tx_id"))
  edata<-collapseBundles(kexp,byWhich,read.cutoff=read.cutoff)
  fit<-lm.fit(design,log2(1+t(edata)))
  d<-DGEList(counts=edata)
  d<-calcNormFactors(d)
  cpm.norm<-cpm(d,normalized.lib.sizes=TRUE,log=FALSE)
  batch.cpm<-removeBatchEffect(log2(1+cpm.norm),batch=batch,design=design)
  be_fit<-lm.fit(design,t(batch.cpm))
  par(mfrow=c(2,2))
  hist(fit$coefficients[2,],col=2,breaks=100,xlab="Fitted without Batch Correction Coefficients",main="Fitted Value Distributions")
  hist(be_fit$coefficients[2,],col=2,breaks=100, xlab="Fitted with Batch Correction Coefficients",main="Fitted Value Distributions")
  plot(fit$coefficients[2,],be_fit$coefficients[2,],col=2,xlab="Without Batch Corrections Fit",ylab="Batch Corrected Fitted Values",main="Batch Correction Fitted Plot")

 par(mfrow=c(2,2))
 pdf(paste0(byLevel,"_",byWhich,"_batchEffectNormalized.pdf"))
   hist(fit$coefficients[2,],col=2,breaks=100,xlab="Fitted without Batch Correction Coefficients",main="Fitted Value Distributions")
  hist(be_fit$coefficients[2,],col=2,breaks=100, xlab="Fitted with Batch Correction Coefficients",main="Fitted Value Distributions")
  plot(fit$coefficients[2,],be_fit$coefficients[2,],col=2,xlab="Without Batch Corrections Fit",ylab="Batch Corrected Fitted Values",main="Batch Correction Fitted Plot")
 dev.off()

  ###correct out NA values
  id<-which(is.na(cpm))
    cpm[id]<-1.0
    rid<-which(is.na(rpm))
    rpm[rid]<-1.0


  if(byLevel=="repeat"){
  batch.cpm<-batch.cpm[!grepl("^ENS",rownames(batch.cpm)),]
  batch.cpm<-batch.cpm[!grepl("^ERCC",rownames(batch.cpm)),]
  }else if(byLevel=="gene"){
  batch.cpm<-batch.cpm[grepl("^ENS",rownames(batch.cpm)),]
 }
 return(batch.cpm)
}
