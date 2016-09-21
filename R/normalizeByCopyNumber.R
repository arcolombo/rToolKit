#' @title each repeat has lots of copy number by nature, so we normalize by copy number to have expression count per copy
#' @description divides the repeat count by copy number
#' @param kexp a kallistoExperiment 
#' @export
#' @return a matrix or data frame of copy normalized counts and tpm
normalizeByCopyNumber<-function(kexp,species=c("Homo.sapiens","Mus.musculus")){
  
  if(species=="Homo.sapiens"){
   data(repeatCopyNumber.hg19,package="repeatToolKit")
   rpm<-collapseBundles(kexp)
   x<-rpm[rownames(rpm)%in%rownames(repeatCopyNumber.hg19),] 
  y<-repeatCopyNumber.hg19[rownames(repeatCopyNumber.hg19)%in%rownames(x),]
 id<-match(rownames(x),rownames(y))
  stopifnot(rownames(x)==rownames(y[id,]))
  copy.norm.cpm<-x/y[id,2]
 metadata(kexp)$copy.norm.cpm<-copy.norm.cpm
   tpm<-collapseTpm(kexp)
   tx<-tpm[rownames(tpm)%in%rownames(repeatCopyNumber.hg19),]
  ty<-repeatCopyNumber.hg19[rownames(repeatCopyNumber.hg19)%in%rownames(tx),]
 tid<-match(rownames(tx),rownames(ty))
  stopifnot(rownames(tx)==rownames(ty[tid,]))
  copy.norm.tpm<-tx/ty[tid,2]
   metadata(kexp)$copy.norm.tpm<-copy.norm.tpm
   } else if(species=="Mus.musclus"){
      data(mouseCopyNumber.mm10,package="repeatToolKit")
   rpm<-collapseBundles(kexp)
   x<-rpm[rownames(rpm)%in%rownames(mouseCopyNumber.mm10),]
  y<-mouseCopyNumber.mm10[rownames(mouseCopyNumber.mm10)%in%rownames(x),]
 id<-match(rownames(x),rownames(y))
  stopifnot(rownames(x)==rownames(y[id,]))
  copy.norm.cpm<-x/y[id,2]
  metadata(kexp)$copy.norm.cpm<-copy.norm.cpm
   tpm<-collapseTpm(kexp)
   tx<-tpm[rownames(tpm)%in%rownames(mouseCopyNumber.mm10),]
  ty<-mouseCopyNumber.mm10[rownames(mouseCopyNumber.mm10)%in%rownames(tx),]
 tid<-match(rownames(tx),rownames(ty))
  stopifnot(rownames(tx)==rownames(ty[tid,]))
  copy.norm.tpm<-tx/ty[tid,2]
   metadata(kexp)$copy.norm.tpm<-copy.norm.tpm
  }

  message("copy normalized repeats is stored under metadata")
 return(kexp)
} ##main
