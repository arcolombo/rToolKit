library(qusage)
library(arkas)
library(repeatToolKit)
library(edgeR)

load("~/Documents/Arkas-Paper-Data/AML-bonemarrow-LSCs/patient-plot-data/wgcna_data/repeat_short_list_biotypes/amlX.RData")

kexp<-kexpByTrio(amlX)
phsc.lsc<-kexp2Group(kexp,comparison="pHSC",control="LSC")

counts<-collapseBundles(phsc.lsc,"gene_name",read.cutoff=6)

 counts<-counts[,which(colSums(counts)>1)]
   dge<-DGEList(counts=counts)
   dge<-calcNormFactors(dge)
   expr<-cpm(dge,normalized.lib.sizes=TRUE,log=FALSE)
  expr2<-log2(1+expr)

cN<-colnames(expr2)
cN<-strsplit(cN,"_")
labels<-unlist(lapply(cN,function(x) x[1]))

contrast<-"pHSC-LSC"

 pairs<-unlist(lapply(cN,function(x) x[2]))
   pairs.id<-match(toupper(pairs),toupper(pairs))
geneSets<-read.gmt("~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/c5.all.v5.1.symbols.gmt")

 qs.pHSC.results<-qusage(expr2,labels,contrast,geneSets,pairVector=pairs.id)
  ##pathway indices : inflam  903, immune 394,410 ,wound 831, mapK 865,494, 
 plot(qs.pHSC.results)
 x<-qsTable(qs.pHSC.results,number=950)
  for(i in c(903,394,410,831,865,494)){
  #in
  par(mfrow=c(1,2),cex.main=0.7)
  plotDensityCurves(qs.pHSC.results,path.index=i,xlim=c(-2,2),col=3,main=paste0(x[which(rownames(x)== i),1]," pHSC-LSC")  ) 
  plotCIsGenes(qs.pHSC.results,path.index=i, main=paste0(x[which(rownames(x)== i),1]," pHSC-LSC") ) 
  readkey()
 ##its own page
 plotGeneSetDistributions(qs.pHSC.results,path.index=i ,groupLabel="pHSC-LSC",main=paste0(x[which(rownames(x)== i),1]," pHSC-LSC")  )
###
readkey()
}
save(qs.pHSC.results,file="pHSC.vs.LSC.qusageResults_CPM.RData",compress=TRUE)
