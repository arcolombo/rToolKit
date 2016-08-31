#' @title stage wise analysis of clonal expression using TMM normalization
#' @description stage level differential expression for global monotonicity analysis using TMM library normalization
#' @param kexp a kallisto expriment at the repeat level stage  
#' @param stage1 first stage
#' @param stage2 second stage
#' @param stage3 third stage
#' @param read.cutoff integer floor bundle count
#' @param byLevel a character of the level for DE analysis
#' @param p.value numeric sig leve
#' @param fold.cutoff floor for logFC 
#' @import limma
#' @import arkas
#' @import beeswarm
#' @export
stageWiseAnalysis<-function(kexp,stage1="pHSC",stage2="LSC",stage3="Blast",read.cutoff=1, byLevel=c("transcript","tx_biotype","gene_biotype"),p.value=0.05,fold.cutoff=0.4 ) {
  ##TMM Normalization
 #must form kexp2Group  for LSC/pHSC
 LSC.v.pHSC<-kexp2Group(kexp,comparison=stage2,control=stage1) 
 rwa.stage1<-repeatWiseAnalysis(LSC.v.pHSC,
              design=metadata(LSC.v.pHSC)$design,
              how="cpm",
              read.cutoff=read.cutoff,
              species="Homo.sapiens",
              adjustBy="BH",
              fold.cutoff=fold.cutoff)
  top.stage1<-topTable(rwa.stage1$fit,n=nrow(LSC.v.pHSC))
  top.stage1.PVal<-top.stage1[which(top.stage1$P.Value<=p.value),]
  top.stage1.adj<-top.stage1[which(top.stage1$adj.P.Val<=p.value),]

 #must form kexp2Group for Blast/LSC
 Blast.v.LSC<-kexp2Group(kexp,comparison=stage3,control=stage2)
 rwa.stage2<-repeatWiseAnalysis(Blast.v.LSC,
              design=metadata(Blast.v.LSC)$design,
              how="cpm",
              read.cutoff=read.cutoff,
              species="Homo.sapiens",
              adjustBy="BH",
              fold.cutoff=fold.cutoff)
  top.stage2<-topTable(rwa.stage2$fit,n=nrow(Blast.v.LSC))
  top.stage2.PVal<-top.stage2[which(top.stage2$P.Value<=p.value),]
  top.stage2.adj<-top.stage2[which(top.stage2$adj.P.Val<=p.value),]

 
 Blast.v.pHSC<-kexp2Group(kexp,comparison=stage3,control=stage1)
 rwa.stage3<-repeatWiseAnalysis(Blast.v.pHSC,
              design=metadata(Blast.v.pHSC)$design,
              how="cpm",
              read.cutoff=read.cutoff,
              species="Homo.sapiens",
              adjustBy="BH",
              fold.cutoff=fold.cutoff)
  top.stage3<-topTable(rwa.stage3$fit,n=nrow(Blast.v.pHSC))
  top.stage3.PVal<-top.stage3[which(top.stage3$P.Value<=p.value),]
  top.stage3.adj<-top.stage3[which(top.stage3$adj.P.Val<=p.value),]
 
   if(byLevel=="transcript"){
  ##based on PVal
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(asinh(top.stage1.PVal$logFC), main=expression(paste(Delta,"(pHSC,LSC)  TMM Norm")),xlab=paste0("LSC-pHSC p.val",p.value) )
  par(new=TRUE)
  boxplot(asinh(top.stage1.PVal$logFC),medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  beeswarm(asinh(top.stage2.PVal$logFC),main=expression(paste(Delta,"(LSC,Blast)  TMM Norm")),xlab=paste0("Blast-LSC p.val",p.value))
  par(new=TRUE)
  boxplot(asinh(top.stage2.PVal$logFC),medcol="red",boxcol="red",whiskcol="red")
  readkey()
 
 ##based on adj.Pval  FDR "BH"
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(asinh(top.stage1.adj$logFC), main=expression(paste(Delta,"(pHSC,LSC)  TMM Norm BH")),xlab=paste0("LSC-pHSC p.val",p.value) )
  par(new=TRUE)
  boxplot(asinh(top.stage1.adj$logFC),medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
 
  beeswarm(asinh(top.stage2.adj$logFC),main=expression(paste(Delta,"(LSC,Blast)  TMM Norm BH")),xlab=paste0("Blast-LSC p.val",p.value))
  par(new=TRUE)
  boxplot(asinh(top.stage2.adj$logFC),medcol="red",boxcol="red",whiskcol="red")
  readkey() 
 
  ###plot the frequency
  plotFrequency(LSC.v.pHSC,topNames=rownames(top.stage1.PVal),topDE=top.stage1.PVal,whichDelta="delta1",p.cutoff=p.value,isAdjusted=FALSE)
  plotFrequency(LSC.v.pHSC,topNames=rownames(top.stage1.adj),topDE=top.stage1.adj,whichDelta="delta1",p.cutoff=p.value,isAdjusted=TRUE)
  ###delta2 PVal and adj.P.Val
 plotFrequency(Blast.v.LSC,topNames=rownames(top.stage2.PVal),topDE=top.stage2.PVal,whichDelta="delta2",p.cutoff=p.value)
  plotFrequency(Blast.v.LSC,topNames=rownames(top.stage2.adj),topDE=top.stage2.adj,whichDelta="delta2",p.cutoff=p.value,isAdjusted=TRUE)
  } else if(byLevel=="tx_biotype") {
  #find DE of Tx Dbiotypes for each stage
     delta1_tx_fit<-fitTxBiotypes(LSC.v.pHSC,design=metadata(LSC.v.pHSC)$design)
    delta1_tx_top<- topTable(delta1_tx_fit$fit,n=nrow(LSC.v.pHSC),p.value=p.value)
     
     delta2_tx_fit<-fitTxBiotypes(Blast.v.LSC,design=metadata(Blast.v.LSC)$design)
     delta2_tx_top<- topTable(delta2_tx_fit$fit,n=nrow(Blast.v.LSC),p.value=p.value)

  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(asinh(delta1_tx_top$logFC), main=expression(paste(Delta,"(pHSC,LSC) DE Tx Biotype")),xlab=paste0("LSC-pHSC p.val",p.value),ylab="TMM Norm BH" )
  par(new=TRUE)
  boxplot(asinh(delta1_tx_top$logFC),medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  beeswarm(asinh(delta2_tx_top$logFC),main=expression(paste(Delta,"(LSC,Blast) DE Tx Biotype")),xlab=paste0("Blast-LSC p.val",p.value),ylab="TMM Norm BH")
  par(new=TRUE)
  boxplot(asinh(delta2_tx_top$logFC),medcol="red",boxcol="red",whiskcol="red")
  readkey()
 } else { 
   #gene biotype DE 
  message("not supporting Gene biotype expression")


 } ##DE lists



  ##FIX ME: Repeat this analysis with CQN and edgeR
 


} ###{ main
