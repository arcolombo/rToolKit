#' @title calculates the connectivity of the Module eigenvalues across genes
#' @description the module eigenvalues are defined as the variability of a given module, and the kME is the connectivity of genes to the Module eigenvalues.  if genes are highly connected to the moduleseigen values then they will behave connectedly to the gene corresponding to the ME.
#' @param datExpr1 expression data 
#' @param ME_G1  module eigenvalues for the expression data
#' @param colorsG1 the color assignment
#' @param modulesG1  the partition of the expression data
#' @param species  Homo.sapiens
#' @param main  gene or Repeat
kmeTable<-function(datExpr1=NULL,ME_G1=NULL,colorsG1=NULL,modulesG1=NULL,species="Homo.sapiens",main="Gene"){
  main<-match.arg(main,c("Gene","Repeat"))
 geneModuleMembership1 = signedKME(t(datExpr1), ME_G1)
  colnames(geneModuleMembership1)=paste("PC",colorsG1,".cor",sep="");
  MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExpr1)[[2]]);
  colnames(MMPvalue1)=paste("PC",colorsG1,".pval",sep="");
  Gene= rownames(datExpr1)
  kMEtable1 = cbind(Gene,Gene,modulesG1)
  for (i in 1:length(colorsG1)){
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
  }
 colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

   return(kMEtable1)
}

