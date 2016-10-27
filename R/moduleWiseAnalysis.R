#' @title this analyzes correlations of repeat and gene modules of fixed interest
#' @description after running wgcna and wrcna we are interested in finding associaitons between one fixed gene with all the repeat modules.  given a fixed module of interest we correlate and test for comparability between every repeat module.  If there is no contrast between the gene modules, i.e. qusage contrasts are not present, then we use the wgcna_goEnrich caller
#' @param kexp this should be the full stage kexp, could even be a tcga kexp
#' @param lnames is the wgcna.R output and rnames is the wrcna.R call
#' @param wgcnaDbName this is the sql db of the wgcna
#' @param wrcnaDbName sql db of the wrcna repeat db
#' @param qusageDbName this is the sql db of a qusage db. not required if there is a missing qusage db one can call goEnrich for preliminary analyses
#' @param geneModules this should be single gene module of interest.  we can examine one gene module's association to all the repeat modules.  so geneModule should be a single color character
#' @param repeatModules this should be the character of all repeat modules
#' @param repeat.cutoff the cpm floor threshold
#' @param how  either cpm or tpm, tpm recommended for heatmap visualizations
#' @param openDevice boolean if true then the images are opened on a device along with printing to pdf.  should be false for exploring all colors
#' @import WGCNA
#' @import edgeR
#' @import pvclust
#' @import ComplexHeatmap
#' @export
moduleWiseAnalysis<-function(kexp,lnames,rnames,wgcnaDbName="wgcnaDbLite.cpm.sqlite",wrcnaDbName="wrcnaDbLite.cpm.sqlite",qusageDbName="qusageDbLite.cpm.sqlite", geneModules=c("darkturquoise","red","grey60","pink","darkred","royalblue","green","magenta","brown"),repeatModules=c("blue","brown","grey","turquoise"),read.cutoff=2,how="cpm",enrichmentCaller=c("go","qusage"),comparison1="phsc",comparison2="blast",openDevice=FALSE ){
  ##FIX ME: generically handle repeat modules of various sizes.
  how<-match.arg(how,c("cpm","tpm"))
  ## gather all ENGIDs for each module.
  enrichmentCaller<-match.arg(enrichmentCaller,c("go","qusage"))
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  rownames(MEs)<-colnames(kexp)
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  stopifnot(lnames[["byWhich"]]=="gene")

  rbwModuleColors<-rnames[["moduleColors"]]
  rMEs<-rnames[["MEs"]]
  rdatExpr<-rnames[["datExpr"]]
  rdatTraits<-rnames[["datTraits"]]
  annot<-rnames[["annot"]]
  stopifnot(rnames[["byWhich"]]=="repeat")
  allcolors<-dbListTables(dbconn(wgcnaDbLite(wgcnaDbName)))
  geneModules<-match.arg(geneModules,allcolors)
  repeatColors<-dbListTables(dbconn(wgcnaDbLite(wrcnaDbName)))
  stopifnot(all(repeatModules%in%repeatColors)==TRUE)
  
  
  gene.color<-modulesBy(wgcnaDbLite(wgcnaDbName),Module.color=geneModules)
  repeat.module.list<-lapply(repeatModules,function(x) modulesBy(wgcnaDbLite(wrcnaDbName),Module.color=x))
  names(repeat.module.list)<-repeatModules

  if(enrichmentCaller=="qusage"){
  ## print enrichment data to screen and write out to CSV
  pHSC<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast=comparison1)
  if(nrow(pHSC)>0){
  pHSC<-data.frame(pHSC,contrast=comparison1,color=geneModules)
  }
  blast<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast=comparison2)  
  if(nrow(blast)>0){
  blast<-data.frame(blast,contrast=comparison2,color=geneModules)
  }
  if(nrow(pHSC)>0 &&nrow(blast)>0){
  pathWays<-rbind(pHSC,blast)
  } else if(nrow(pHSC)==0 &&nrow(blast)>0){
   pathWays<-blast
  } else if(nrow(pHSC)>0 && nrow(blast)==0){
   pathWays<-pHSC
  } else{
  pathWays<-data.frame(pathway=geneModules,data="NA")
  }
  print(pathWays)
  write.csv(pathWays ,file=paste0(geneModules,"_",how,"_moduleWiseAnalysis.csv"))
  readkey()
  }else{
  #if no qusage contrasts, call goEnrich and print out the fixed color module
  goenriched<-goEnrich(wgcnaDbLite(wgcnaDbName),Module.color=geneModules)
  print(goenriched)
  write.csv(as.data.frame(goenriched) ,file=paste0(geneModules,"_",how,"_GOmoduleWiseAnalysis.csv"))
  readkey()
  } 
  ##find kexp expression values for each
   ##
 if(how=="cpm"){
  cpm<-collapseBundles(kexp,"gene_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  cpm<-cpm[grepl("^ENS",rownames(cpm)),]
  rexp<-findRepeats(kexp)
  rpm<-collapseBundles(rexp,"tx_id",read.cutoff=read.cutoff)
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
   d<-DGEList(counts=cpm)
  cpm.norm<-cpm(d,normalized.lib.sizes=TRUE)
  cpm<-cpm.norm
  cpm.norm<-NULL
  rd<-DGEList(counts=rpm)
  rdm.norm<-cpm(rd,normalized.lib.sizes=TRUE)
  rpm<-rdm.norm
  rdm.norm<-NULL
 } else{
  cpm<-collapseTpm(kexp,"gene_id",read.cutoff=read.cutoff)
  cpm<-cpm[!grepl("^ERCC",rownames(cpm)),]
  cpm<-cpm[grepl("^ENS",rownames(cpm)),]
  rexp<-findRepeats(kexp)
  rpm<-collapseTpm(rexp,"tx_id",read.cutoff=read.cutoff)
  rpm<-rpm[!grepl("^ERCC",rownames(rpm)),]
 }
  #####

  ###now subset geneExpression for module (gene/repeat) level
  
  cpm.id<-match(gene.color$gene_id,rownames(cpm))
  module.cpm<-cpm[cpm.id,] 
  module.cpm<-module.cpm[!is.na(rownames(module.cpm)),]
  stopifnot(all(rownames(module.cpm)%in%gene.color$gene_id))

  color.ID<-which(paste0("ME",geneModules)==colnames(MEs))
 ### barplot of sample module eigen values for a given module
  dev.new()
  par(mar=c(4,10,4,4))
  barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.8,las=1,main=paste0(geneModules," gene Module Eigen"),space=1)
 
  for(colR in repeatModules){
  dev.new()
    par(mar=c(4,10,4,4))
   color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  ME.id<-grep(paste0("ME",colR),colnames(rMEs))
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
  barplot(rMEs[,ME.id],horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0(colR," Repeat Module EigenValue cor=", bicor(rMEs[,ME.id],MEs[,color.ID])),space=1)
   }
  readkey()
  ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
  if(openDevice==TRUE){
  ###heatmap of module genes
 dev.new()
  print(Heatmap(asinh(module.cpm),column_title=paste0(geneModules," Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=5) ) )
 ##heatmap of repeat modules
  dev.new()
   for(colR in repeatModules){
  color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
  dev.new()
  print(Heatmap(asinh(module.color),column_title=paste0(colR," Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.color)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.color)],name="gene_biotype"))
  }
 }
 

 
 ### PRINT TO PDF########3
  pdf(paste0(geneModules,"_",how,"_",enrichmentCaller,"_EigenValueAnalysis.pdf"))
  par(mar=c(4,10,4,4))
   barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.6,las=1,main=paste0(geneModules," gene Module Eigen"),space=1)
    for(colR in repeatModules){
    par(mar=c(4,10,4,4))
   color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  ME.id<-grep(paste0("ME",colR),colnames(rMEs))
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
  barplot(rMEs[,ME.id],horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0(colR," Repeat Module EigenValue cor=", bicor(rMEs[,ME.id],MEs[,color.ID])),space=1)
   }
 ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
print( Heatmap(asinh(module.cpm),column_title=paste0(geneModules," Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6)))
 ##heatmap of repeat modules
   for(colR in repeatModules){
  color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
  print(Heatmap(asinh(module.color),column_title=paste0(colR," Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.color)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.color)],name="gene_biotype"))
  }
 dev.off()

} ##main
