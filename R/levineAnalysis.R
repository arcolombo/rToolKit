#' @title analyzes levine data
#' @description analyzes the levine kexp by GMP Mx-Cre+ (interferon activated) (3) vs Mx-Cre (neg) Cntrl.  
#' @import arkas
#' @export
#' @return a kexp analysis
levineAnalysis<-function(kexp){

  gmp<-kexp[,grep("GMP",colnames(kexp))]
  design<-matrix(nrow=ncol(gmp),ncol=2)
  rownames(design)<-colnames(gmp)
  design[,1]<-1
  design[grep("pos",rownames(design)),2]<-1
  design[grep("neg",rownames(design)),2]<-0
  colnames(design)<-c("Intercept","Pos.v.Neg")
  design<-design[order(design[,2]),]
  metadata(gmp)$design<-design
  lsk<-kexp[,grep("LSK",colnames(kexp))]
  design2<-matrix(nrow=ncol(lsk),ncol=2)
  rownames(design2)<-colnames(lsk)
  design2[,1]<-1
  design2[grep("pos",rownames(design2)),2]<-1
  design2[grep("neg",rownames(design2)),2]<-0
  colnames(design2)<-c("Intercept","Pos.v.Neg")
  design2<-design2[order(design2[,2]),]
  metadata(lsk)$design<-design2

###gmp analysis
   gmp<-annotateFeatures(gmp,"transcript")
   lsk<-annotateFeatures(lsk,"transcript")
  gmp_meta<-kexpAnalysis(gmp,byWhich="repeat",adjustBy="none",read.cutoff=3)
  boxplot(gmp_meta$top$logFC,main="Overall LogFC Mx pos .v Mx Neg GMP")

  lsk_meta<-kexpAnalysis(lsk,byWhich="repeat",adjustBy="BH",read.cutoff=3)
  boxplot(lsk_meta$top$logFC,main="Overall LogFC Mx pos .v Mx Neg LSK")




##LSK analysis



}
