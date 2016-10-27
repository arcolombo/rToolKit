#' @title each repeat has lots of copy number by nature, so we normalize by copy number to have expression count per copy
#' @description divides the repeat count by copy number
#' @param kexp a kallistoExperiment 
#' @import arkas
#' @export
#' @return a matrix or data frame of copy normalized counts and tpm
normalizeByCopyNumber<-function(kexp,species=c("Homo.sapiens","Mus.musculus")){
 kexp<-findRepeats(kexp) 
  if(species=="Homo.sapiens"){
   data(repeatCopyNumber.hg19,package="repeatToolKit")
   rpm<-collapseBundles(kexp)
   x<-rpm[rownames(rpm)%in%rownames(repeatCopyNumber.hg19),] 
  stopifnot(nrow(x)>0)
  y<-repeatCopyNumber.hg19[rownames(repeatCopyNumber.hg19)%in%rownames(x),]
 id<-match(rownames(x),rownames(y))
  stopifnot(rownames(x)==rownames(y[id,]))
  copy.norm.cpm<-x/y[id,2]
 metadata(kexp)$copy.norm.cpm<-copy.norm.cpm

  ###for TPM we need the est_counts/copyNumber and eff length
   effLength<-eff_length(kexp)
   idEf<-match(rownames(x),rownames(effLength))
   efflng<-effLength[idEf,]
   stopifnot(rownames(efflng)==rownames(copy.norm.cpm))
   tpm.cn<-calcTpm(copy.norm.cpm,efflng)   
   metadata(kexp)$copy.norm.tpm<-tpm.cn
   } else if(species=="Mus.musclus"){
      data(mouseCopyNumber.mm10,package="repeatToolKit")
   rpm<-collapseBundles(kexp)
   x<-rpm[rownames(rpm)%in%rownames(mouseCopyNumber.mm10),]
  y<-mouseCopyNumber.mm10[rownames(mouseCopyNumber.mm10)%in%rownames(x),]
 id<-match(rownames(x),rownames(y))
  stopifnot(rownames(x)==rownames(y[id,]))
  copy.norm.cpm<-x/y[id,2]
  metadata(kexp)$copy.norm.cpm<-copy.norm.cpm

  ##FOR TPM you should not merely divide by copyNumber, but must remap to TPM in the end####   
   effLength<-eff_length(kexp)
   idEf<-match(rownames(x),rownames(effLength))
   efflng<-effLength[idEf,]
   stopifnot(rownames(efflng)==rownames(copy.norm.cpm))
   tpm.cn<-calcTpm(copy.norm.cpm,efflng)     
   metadata(kexp)$copy.norm.tpm<-tpm.cn
 }

  message("copy normalized repeats is stored under metadata")
 return(kexp)
} ##main
