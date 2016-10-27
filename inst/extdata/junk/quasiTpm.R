#' @title quasi normalization of TPM by bundling
#' @description resulting quasi-normalization of TPM of a bundle numbers should now have variances which behave similarly to those of counts. this may be good enough to use the bundle wise TPM quanties for batch normalization and other shenanigans, then multiply through by the original proportions in each sample to retrieve the corrected per-transcript values to do unsupervised work.  this might also address aspects of single cell noise.
#' @param kexp a stage level kexp
#' @param ensg an EnsgID
#' @param sample a sample column
#' @export
#' @return normalized tpm
quasiTpm<-function(kexp,ensg=ensgID,sample=sampleID){
  
  sample_bundleDF<-.findBundle(kexp,ensg=ensgID,sample=sampleID)
  ##first must estimate s_i  
  total_mass<-sum( sample_bundleDF$est_counts/sample_bundleDF$eff_len)
  L_g<-.medEffectiveLength(kexp,ensg=ensgID)
  ##S is the set of all expression features in a bundle  
  s_i<-sum(.scaleSampleFactor(kexp,ensg=ensgID,sample=sampleID))
    ###FIX ME:: s_i needs to be summed across txs , s_i returns a vector for each tanscript in bundle,  need aggregated data?
   sample_quasi_tpm<- (L_g/s_i)/total_mass

} ##main

.findBundle<-function(kexp,ensg=NULL,sample=NULL){
stopifnot(is.null(ensg)==FALSE)
txs<-rowRanges(kexp)[which(rowRanges(kexp)$gene_id==as.character(ensg))]
tx_effective_length<-assays(kexp)$eff_length[txs$tx_id,sample]
sample_name<-colnames(kexp)[sample]
estCounts<-counts(kexp)[txs$tx_id,sample]
tpm<-tpm(kexp)[txs$tx_id,sample]
 df<-data.frame(txs=txs$tx_id,
                eff_len=tx_effective_length,
                est_counts=estCounts,
                ensg=txs$gene_id,
                tpm=tpm,
                sample=sample_name,
                stringsAsFactors=FALSE)

return(df)
}


.medEffectiveLength<-function(kexp,ensg=NULL ) {
txs<-rowRanges(kexp)[which(rowRanges(kexp)$gene_id==as.character(ensg))]
tx_effective_length<-assays(kexp)$eff_length[txs$tx_id,]
return(median(tx_effective_length))
}

.scaleSampleFactor<-function(kexp,ensg=NULL,sample=NULL){
txs<-rowRanges(kexp)[which(rowRanges(kexp)$gene_id==as.character(ensg))]
sample_total_count<-sum(counts(kexp)[txs$tx_id,sampleID])
tx_count<-counts(kexp)[txs$tx_id,sampleID]

 s_i<-tx_count/sample_total_count
 return(s_i)
}
