#' @title this will plot a cluster of eigen values and pca of module eigenvalues
#' @description by clustering moduleEigenvalues one can examine how close other modules are to one another.  if to module eigenvalues are branched in the same tree, then their functions pathway are likely to be similiar.  Further modules with similiar eigenvalues can even be merged.
#' @param lnames this is the wgcna.cpm.RData returned from wgcna.R
#' @import WGCNA
#' @export
modulePCA<-function(lnames){

  ##FIX ME:  it may be a solution to cluster repeatModule eigenvalues with gene module eigenvalues to see which repeat modules fall closest to gene modules?  this may help with the hub-connectivity analysis?

  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  how<-lnames[["how"]]
  byWhich<-lnames[["byWhich"]]
 annot<-lnames[["annot"]]
 stopifnot(lnames[["byWhich"]]=="gene")
  ###pca plot of Module Gene eigen values, and repeat Eigenvalues
  ME_1A<-MEs
  distPC1A<-1-abs(cor(ME_1A,use="p"))
  distPC1A<-ifelse(is.na(distPC1A),0,distPC1A)
  pcTree1A<-hclust(as.dist(distPC1A),method="a")
  MDS_1A<-cmdscale(as.dist(distPC1A),2)
  colorsA1<-names(table(bwModuleColors))
  dev.new()
  par(mfrow=c(1,2), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="ME Hier.Cluster",sub="")
plot(MDS_1A, col= colorsA1,main="MDS plot", cex=2, pch=19)
 readkey()
  pdf(paste0(byWhich,"_",how,"_ClusterMDS.pdf"),width=12,height=7)
 par(mfrow=c(1,2), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="ME Hier.Cluster",sub="")
plot(MDS_1A, col= colorsA1,main="MDS plot", cex=2, pch=19)
 dev.off()

cat("done.\n")
} #main
