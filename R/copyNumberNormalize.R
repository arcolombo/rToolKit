#' @title normalizes the copy number for repeat biotypes family expression
#' @description averages the family biotype expression and normalizes with the average copy number
#' @import GenomicRanges
#' @return a matrix of copy normalized repeat counts
#' @export
copyNumberNormalize<-function(kexp,bundleID="tx_biotype",read.cutoff=1,discardjoined=TRUE){

  message("For the time being, only summing of bundles is supported")

  bundleable <- !is.na(mcols(rowRanges(kexp))[[bundleID]])
  feats <- rowRanges(kexp)[bundleable]
  cts <- split.data.frame(counts(kexp)[bundleable, ], mcols(feats)[[bundleID]])
  cts <- cts[ sapply(cts, function(x) max(x) >= read.cutoff) ]
  cts <- lapply(cts, function(x) x[ rowSums(x) >= read.cutoff, ])
  bundled <- do.call(rbind, lapply(cts,
                                   function(x)
                                     if(!is.null(nrow(x))) colSums(x) else x))


 ###bundling copy number

 copies <- !is.na(mcols(rowRanges(kexp))[["copyNumber"]])
  copy.feats <- rowRanges(kexp)[copies]
  CN<-as.data.frame(rowRanges(kexp)$copyNumber[copies])
  rownames(CN)<-rowRanges(kexp)$tx_id[copies]
  
  copy.cts <- split.data.frame(CN, mcols(feats)[[bundleID]])
  copy.cts <- copy.cts[ sapply(copy.cts, function(x) max(x) >= read.cutoff) ]
  copy.total <- lapply(copy.cts, function(x)nrow(x))
  copy.bundled <- do.call(rbind, lapply(copy.cts,
                                   function(x)
                                     if(!is.null(nrow(x))) colSums(x) else x))

  id<-match(rownames(bundled),rownames(copy.bundled))
  final<-copy.bundled[id,]
  outputs<-bundled/final
###
  if (discardjoined) {
    return(outputs[ grep(";", invert=TRUE, rownames(outputs)), ])
  } else {
    return(outputs)
  }


}
