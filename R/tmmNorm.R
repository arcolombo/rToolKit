#' @title TMM Normalization KallistoExperiment
#' @description input a KallistoExperiment and TMM the library size to normalize sequence depth and store the normalized values in the kexp container. The library normalized assays are added to the kallistoExperiment object which can be used in a patientPlot
#' @param kexp a kallistoExperiment
#' @importFrom edgeR DGEList calcNormFactors cpm rpkm
#' @export
#' @return a kallistoExperiment 
tmmNorm<-function(kexp){

  expr<-DGEList(counts=counts(kexp))
  expr<-calcNormFactors(expr)
  print(expr$samples)
  expr_cpm<-cpm(expr,log=FALSE)
  stopifnot(dim(expr_cpm)==dim(kexp))
  assays(kexp)$cpm_library_normalized<-expr_cpm
 
  stopifnot(all(rownames(expr$counts)==rownames(eff_length(kexp)))==TRUE)
  expr_rpkm<-rpkm(expr,log=FALSE,gene.length=eff_length(kexp))
  stopifnot(dim(expr_rpkm)==dim(kexp))
  assays(kexp)$rpkm_library_normalized<-expr_rpkm
  return(kexp) 
}
