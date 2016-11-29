#' @title tests the connectivity in gene pathways
#' @description tests the significance of connectivity of genes in a given pathway resulting in significance tests of connected genes of pathway.
#' @import WGCNA

connectivityWiseAnalysis<-function(lnames,MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"), pathway=NULL, dbname="wgcnaDBLite.sqlite",qusageDbName="qusageDbLite",gmt.path="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",moduleColor="violet" ){

 bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
 nGenes= ncol(datExpr);
  nGenes= ncol(datExpr);
  nSamples = nrow(datExpr);
  moduleTraitCor = bicor(MEs, datTraits, use = "all.obs");
   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

geneModuleMembership1<-signedKME(datExpr,MEs)
  geneSets<-read.gmt(paste0(gmt.path,"/",geneSet))

###subset the kME set by the pathway of interest
##grab the enriched pathway

##the pathway of interest is in gene_name so maybe need to collapse kME by gene name

###need to collapse kME by gene name?

for(i in c(1,2,3,4,5,6,8,9,10,12,14,16,17,21,23,24,25,27,31,32,33,39)){+ write.csv(sort(unique(feats[feats$gene_id%in%rownames(kme1)[which(abs(kme1[,i])>=0.80)]]$gene_name)),file=paste0(rownames(moduleTraitCor)[i],"_",i,"_topKME_.80.csv"))}



}
