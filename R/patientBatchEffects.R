#' @title this analyzes the batch effects from specific patients, where each patient ID is considered a batch
#' @description performs a batch correction where each patient is a batch so to remove out individual patient effects to make each donor uniform.  it is critical to put the entire 67 samples kexp into the input.
#' @import sva
#' @import edgeR
#' @import limma
#' @param kexp a full kexp of every sample, not just a stage kexp
#' @param byWhat how to collapse 
#' @param cqnLevel  character either stage or repeats
#' @export
#' @return returns images and a patient batch corrected data set
patientBatchEffects<-function(kexp,byWhat=c("counts","tpm"), cqnLevel=c("stage","repeats"),stage1="pHSC",stage2="LSC",stage3="Blast",p.value=0.05){
  message("the kexp should not be stage level and not just be repeats , use the entire kexp")
 byWhat<-match.arg(byWhat,c("counts","tpm"))
 batchvec<-matrix(nrow=ncol(aml),ncol=1)
 colnames(batchvec)<-"condition.batch"
 rownames(batchvec)<-colnames(kexp)
 ##split by donor ID
 
 rowRanges(kexp)[which(rowRanges(kexp)$biotype_class=="repeat")]$gene_id<-rowRanges(kexp)[ which(rowRanges(kexp)$biotype_class=="repeat")]$tx_id 

 if(byWhat=="counts"){
  inputs<-collapseBundles(kexp,"gene_id",read.cutoff=1)
 } else {
  inputs<-collapseTpm(kexp,"gene_id")
 }


 ##patient specific batch effects
 batchvec<-matrix(nrow=ncol(kexp),ncol=1)
 colnames(batchvec)<-"patient.batch"
 rownames(batchvec)<-colnames(kexp)
 ##split by donor ID
 phenos<-pData(kexp)$ID
 batch<-factor(sapply(strsplit(phenos,"_"),function(x) x[2]))
 for(i in 1:length(batch)){
  batchvec[grep(as.character(batch[i]),rownames(batchvec)),]<-as.character(batch[i])
 }
  pData(kexp)$patient.batch<-as.factor(batchvec)
  psva_patient<-psva(asinh(inputs),pData(kexp)$patient.batch)
  psva_patient_repeats<-psva_patient[!grepl("^ENSG",rownames(psva_patient)),]
  psva_patient_repeats<-psva_patient_repeats[!grepl("^ERCC",rownames(psva_patient_repeats)),]
  colnames(psva_patient_repeats)<-colnames(inputs)
  stageKexp<-kexpByStage(kexp)
  psva_stage<-psva_patient_repeats[,colnames(psva_patient_repeats)%in%colnames(stageKexp)]
  psva_stage<-psva_stage^2
 te2<-HeatmapAnnotation(as.data.frame(pData(kexp)$patient.batch))
 patient_repeat_SVA<-Heatmap(psva_stage,top_annotation=te2,
             name="cpm")
 patient_repeat_SVA2<-Heatmap(rowRanges(kexp)[rownames(psva_stage)]$tx_biotype,
              name="tx_biotype",
              width=unit(5,"mm"))
 patient_repeat_SVA3<-Heatmap(rowRanges(kexp)[rownames(psva_stage)]$gene_biotype,
              name="gene_biotype",
              width=unit(5,"mm"))
 draw(patient_repeat_SVA+patient_repeat_SVA2+patient_repeat_SVA3)
 title(sub="Full Repeat Patient Batch Crctn")
 readkey()

###heatmap expression byMad 
 heatmapByMad(psva_stage,kexp,topAnnoFactors=pData(kexp)$patient.batch,selectK=100,byWhat=byWhat)
 readkey()



###FIX ME: add the treatment batch 
  res<-list()
  reps<-findRepeats(kexp)
  repStage<-kexpByStage(reps)
  LSC.v.pHSC<-kexp2Group(repStage,comparison="LSC",control="pHSC")
  design<-metadata(LSC.v.pHSC)$design
  psva_patient_repeats_cnts<-psva_patient_repeats^2
  dge<-DGEList(counts=psva_patient_repeats_cnts[,colnames(psva_patient_repeats_cnts)%in%rownames(design)]  )
  dge<-calcNormFactors(dge)
  res$design<-design
  res$voomed<-voom(dge,res$design)
  res$fit<-eBayes(lmFit(res$voomed,res$design))
  topTags.none<-topTable(res$fit,adjust.method="none",p=0.05,n=nrow(psva_patient_repeats_cnts),coef=2)
   topTags.bh<-topTable(res$fit,adjust.method="BH",p=0.05,n=nrow(psva_patient_repeats_cnts),coef=2)
#### Group Biotype composition
   plotBatchFrequency(dge$counts,topNames=rownames(topTags.none),topDE=topTags.none,whichDelta="delta1",isAdjusted=FALSE,kexp=repStage)
   plotBatchFrequency(dge$counts,topNames=rownames(topTags.bh),topDE=topTags.bh,whichDelta="delta1",isAdjusted=TRUE,kexp=repStage)

drawBatchHeatmap(repStage,tags=topTags.bh,batchData=dge$counts)

##delta2
 res<-list()
 
  Blast.v.LSC<-kexp2Group(repStage,comparison="Blast",control="LSC")
  design<-metadata(Blast.v.LSC)$design
  dge2<-DGEList(counts=psva_patient_repeats_cnts[,colnames(psva_patient_repeats_cnts)%in%rownames(design)]  )
  dge2<-calcNormFactors(dge2)
  res$design<-design
  res$voomed<-voom(dge2,res$design)
  res$fit<-eBayes(lmFit(res$voomed,res$design))
  topTags.none2<-topTable(res$fit,adjust.method="none",p=0.05,n=nrow(psva_patient_repeats_cnts),coef=2)
   topTags.bh2<-topTable(res$fit,adjust.method="BH",p=0.05,n=nrow(psva_patient_repeats_cnts),coef=2)
#### Group Biotype composition
   plotBatchFrequency(dge2$counts,topNames=rownames(topTags.none2),topDE=topTags.none2,whichDelta="delta2",isAdjusted=FALSE,kexp=repStage)
   plotBatchFrequency(dge2$counts,topNames=rownames(topTags.bh2),topDE=topTags.bh2,whichDelta="delta2",isAdjusted=TRUE,kexp=repStage)

drawBatchHeatmap(repStage,tags=topTags.bh2,batchData=dge2$counts)


### Stage Beeswarm with boxplot
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(asinh(topTags.none$logFC), main=expression(paste(Delta,"(pHSC,LSC) DE Tx Biotype")),xlab=paste0("LSC-pHSC p.val",p.value),ylab="SVA Norm TMM" )
  par(new=TRUE)
  boxplot(asinh(topTags.none$logFC),medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  beeswarm(asinh(topTags.none2$logFC),main=expression(paste(Delta,"(LSC,Blast) DE Tx Biotype")),xlab=paste0("Blast-LSC p.val",p.value),ylab="SVA Norm TMM")
  par(new=TRUE)
  boxplot(asinh(topTags.none2$logFC),medcol="red",boxcol="red",whiskcol="red")
  readkey()
 
  plot.new()
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(asinh(topTags.bh$logFC), main=expression(paste(Delta,"(pHSC,LSC) DE Tx Biotype")),xlab=paste0("LSC-pHSC adj.p.val",p.value),ylab="SVA Norm TMM BH" )
  par(new=TRUE)
  boxplot(asinh(topTags.bh$logFC),medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  beeswarm(asinh(topTags.bh2$logFC),main=expression(paste(Delta,"(LSC,Blast) DE Tx Biotype")),xlab=paste0("Blast-LSC p.val",p.value),ylab="SVA Norm TMM BH")
  par(new=TRUE)
  boxplot(asinh(topTags.bh2$logFC),medcol="red",boxcol="red",whiskcol="red")
  readkey()


}###main
