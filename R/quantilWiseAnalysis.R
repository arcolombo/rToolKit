#' @title splits quartiles of repeat expression
#' @description this will collapse repeats and split by quartile, and plot the patient annotations across quartiles
#' @import gplots
#' @import edgeR
#' @import qusage
quantileWiseAnalysis<-function(tcga,tx.Biotype="Alu",how=c("cpm","tpm")){
 
  how<-match.arg(how,c("cpm","tpm"))
  alu.id<-which(rowRanges(tcga)$tx_biotype==tx.Biotype)
  alus<-tcga[alu.id,]
  alu.tpm<-collapseTpm(alus,bundleID="tx_biotype")
  stats<-summary(alu.tpm)
  q1<-stats[2]
  q2<-stats[3]
  q3<-stats[5] 
  Q1<-sort(alu.tpm[which(alu.tpm<q1)])
  Q1.id<-match(names(Q1),names(alu.tpm))
  q2.min.id<-which(alu.tpm>=q1)
  q2.max.id<-which(alu.tpm<q2)
  Q2<-sort(alu.tpm[intersect(q2.min.id,q2.max.id)])
   q3.min.id<-which(alu.tpm>=q2)
   q3.max.id<-which(alu.tpm<q3)

  Q3<-sort(alu.tpm[intersect(q3.min.id,q3.max.id)]) 
  Q4<-sort(alu.tpm[ which(alu.tpm>=q3) ])
  Q<-c(Q1,Q2,Q3,Q4)
  plot(Q,type='h')
  title(tx.Biotype)
  abline(v=c(length(Q1),length(Q1)+length(Q2),length(Q1)+length(Q2)+length(Q3)))
  readkey()
  mutant.info<-metadata(tcga)$quantileMutants
  q.id<-match(names(Q),colnames(mutant.info))
  mutant.matrix<-mutant.info[,q.id]
  m2<-mutant.matrix
 m2[]<-c(0)[match(m2,"")]
  m2[is.na(m2)]<-1
  
  p53<-grep("p53",rownames(m2),ignore.case=T)
  idh<-grep("idh",rownames(m2),ignore.case=T)
  flt3<-grep("flt3",rownames(m2),ignore.case=T)
  npm1<-grep("npm1",rownames(m2),ignore.case=T)
  runx<-grep("runx",rownames(m2),ignore.case=T)
  dnmt3<-grep("dnmt3",rownames(m2),ignore.case=T)
  tet2<-grep("tet2",rownames(m2),ignore.case=T)
  cebp<-grep("cebp",rownames(m2),ignore.case=T)
  asxl<-grep("askl",rownames(m2),ignore.case=T)
  m.bin<-m2
  class(m.bin)<-"numeric"
  heatmap.2(as.matrix(m.bin[c(p53,idh,flt3,npm1,runx,dnmt3,tet2,cebp,asxl),]),
  dendrogram='none',
   Rowv=FALSE,Colv=FALSE,trace='none',main=tx.Biotype)
  readkey()
 
##split out Q4, and Q1 ####
  quartile<-tcga[,colnames(tcga)%in%c(names(Q1),names(Q4))]
  
if(how=="cpm"){
   quartile.counts<-collapseBundles(quartile,"gene_name")
   dge<-DGEList(counts=quartile.counts)
   dge<-calcNormFactors(dge)
   quartile.expr<-cpm(dge,log=FALSE)
   quartile.counts<-log2(+quartile.expr) ##enirhcment on log2 is required
    } else {
  quartile.counts<-collapseTpm(quartile,"gene_name")
  quartile.counts<-log2(1+quartile.counts)
   }
 
  quartile1.id<-match(names(Q1),colnames(quartile.counts))
  quartile2.id<-match(names(Q4),colnames(quartile.counts))
  quart1<-quartile.counts[,quartile1.id]
  quart2<-quartile.counts[,quartile2.id]
  colnames(quart1)<-paste0("LowerQ_",colnames(quart1))
 colnames(quart2)<-paste0("UpperQ_",colnames(quart2))
  quartF<-cbind(quart1,quart2)
   qusage_run1<-qusageRun(cnts_mt=quartF,MsigDB="h.all.v5.1.symbols.gmt",comparison="UpperQ",control="LowerQ",module=tx.Biotype,paired=FALSE)


  pdf(paste0(tx.Biotype,"mutantplot.pdf"),width=12,height=12)
  heatmap.2(as.matrix(m.bin[c(p53,idh,flt3,npm1,runx,dnmt3,tet2,cebp,asxl),]),
  dendrogram='none',
   Rowv=FALSE,Colv=FALSE,trace='none',main=tx.Biotype)
  plot(Q,type='h')
  title(tx.Biotype)
  abline(v=c(length(Q1),length(Q1)+length(Q2),length(Q1)+length(Q2)+length(Q3)))
  dev.off()
 write.csv(qusage_run1,file=paste0(tx.Biotype,"_qusageTCGA.csv"))
 return(qusage_run1)
 
} #main
