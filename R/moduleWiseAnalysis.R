#' @title this analyzes correlations of repeat and gene modules of fixed interest
#' @description after running wgcna and wrcna we are interested in finding associaitons between gene and repeat modules.  given modules of interest we correlate and test for comparability between each module.  we also analyze the associations of connectivity.
#' @import WGCNA
#' @import edgeR
moduleWiseAnalysis<-function(kexp,lnames,rnames,wgcnaDbName="wgcnaDbLite.cpm.sqlite",wrcnaDbName="wrcnaDbLite.cpm.sqlite",qusageDbName="qusageDbLite.cpm.sqlite", geneModules=c("darkturquoise","red","grey60","pink","darkred","royalblue","green","magenta","brown"),repeatModules=c("blue","brown","grey","turquoise"),read.cutoff=2,how=c("cpm","tpm") ){
  ## gather all ENGIDs for each module.
  
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
  
  

  gene.color<-modulesBy(wgcnaDbLite(wgcnaDbName),Module.color=geneModules)
  repeat.blue<-modulesBy(wgcnaDbLite(wrcnaDbName),Module.color="blue")
  repeat.brown<-modulesBy(wgcnaDbLite(wrcnaDbName),Module.color="brown")
  repeat.grey<-modulesBy(wgcnaDbLite(wrcnaDbName),Module.color="grey")
  repeat.turquoise<-modulesBy(wgcnaDbLite(wrcnaDbName),Module.color="turquoise")

  ## print enrichment data to screen and write out to CSV
  pHSC<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast="phsc")
  if(nrow(pHSC)>0){
  pHSC<-data.frame(pHSC,contrast="phsc",color=geneModules)
  }
  blast<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast="blast")  
  if(nrow(blast)>0){
  blast<-data.frame(blast,contrast="blast",color=geneModules)
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
  ##find kexp expression values for each
   ##FIX ME: add TPM instead?
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

  ##blue
  blue.id<-match(repeat.blue$gene_id,rownames(rpm))
  module.blue<-rpm[blue.id,]
  stopifnot(all(rownames(module.blue)==repeat.blue$gene_id)==TRUE)
  ##brown
  brown.id<-match(repeat.brown$gene_id,rownames(rpm))
  module.brown<-rpm[brown.id,]
  stopifnot(all(rownames(module.brown)==repeat.brown$gene_id)==TRUE)
  ##grey
  grey.id<-match(repeat.grey$gene_id,rownames(rpm))
  module.grey<-rpm[grey.id,]
  stopifnot(all(rownames(module.grey)==repeat.grey$gene_id)==TRUE)
 ##turquiose
  turquoise.id<-match(repeat.turquoise$gene_id,rownames(rpm))
  module.turquoise<-rpm[turquoise.id,]
  stopifnot(all(rownames(module.turquoise)==repeat.turquoise$gene_id)==TRUE)

 color.ID<-which(paste0("ME",geneModules)==colnames(MEs))
 ### barplot of sample module eigen values for a given module
  dev.new()
  par(mfrow=c(3,2),mar=c(4,10,4,4))
  barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.8,las=1,main=paste0(geneModules," gene Module Eigen"),space=1)
 barplot(rMEs$MEblue,horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0("Blue Repeat Module EigenValue cor=", bicor(rMEs$MEblue,MEs[,color.ID])),space=1)
 barplot(rMEs$MEturquoise,horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0("Turquoise Repeat Module EigenValue cor=",bicor(rMEs$MEturquoise,MEs[,color.ID])),space=1)
 barplot(rMEs$MEbrown,horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0("Brown Repeat Module EigenValue cor=",bicor(rMEs$MEbrown,MEs[,color.ID])),space=1)
  barplot(rMEs$MEgrey,horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0("Grey Repeat Module EigenValue cor=",bicor(rMEs$MEgrey,MEs[,color.ID])),space=1)
   readkey()
 
  ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 

  ###heatmap of module genes
 dev.new()
  print(Heatmap(asinh(module.cpm),column_title=paste0(geneModules," Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=5) ) )
 ##heatmap of repeat modules
  dev.new()
  print(Heatmap(asinh(module.blue),column_title=paste0("Blue Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.blue)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.blue)],name="gene_biotype"))
  dev.new()
  print(Heatmap(asinh(module.brown),column_title=paste0("Brown Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.brown)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.brown)],name="gene_biotype"))
 dev.new()
  print(Heatmap(asinh(module.grey),column_title=paste0("Grey Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.grey)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.grey)],name="gene_biotype"))
 dev.new()
  print(Heatmap(asinh(module.turquoise),column_title=paste0("Turquoise Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.turquoise)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.turquoise)],name="gene_biotype"))
 

 
 ### PRINT TO PDF########3
  pdf(paste0(geneModules,"_",how,"_EigenValueAnalysis.pdf"))
  par(mfrow=c(3,2),mar=c(4,10,4,4))
   barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.6,las=1,main=paste0(geneModules," gene Module Eigen"),space=1)
 barplot(rMEs$MEblue,horiz=T,names.arg=colnames(kexp),cex.names=0.6,las=1,main=paste0("Blue Repeat Module EigenValue cor=", bicor(rMEs$MEblue,MEs[,color.ID])),space=1)
 barplot(rMEs$MEturquoise,horiz=T,names.arg=colnames(kexp),cex.names=0.6,las=1,main=paste0("Turquoise Repeat Module EigenValue cor=",bicor(rMEs$MEturquoise,MEs[,color.ID])),space=1)
 barplot(rMEs$MEbrown,horiz=T,names.arg=colnames(kexp),cex.names=0.6,las=1,main=paste0("Brown Repeat Module EigenValue cor=",bicor(rMEs$MEbrown,MEs[,color.ID])),space=1)
  barplot(rMEs$MEgrey,horiz=T,names.arg=colnames(kexp),cex.names=0.6,las=1,main=paste0("Grey Repeat Module EigenValue cor=",bicor(rMEs$MEgrey,MEs[,color.ID])),space=1)
 ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
print( Heatmap(asinh(module.cpm),column_title=paste0(geneModules," Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6)))
 ##heatmap of repeat modules
  print(Heatmap(asinh(module.blue),column_title=paste0("Blue Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.blue)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.blue)],name="gene_biotype"))
 print( Heatmap(asinh(module.brown),column_title=paste0(" Brown Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.brown)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.brown)],name="gene_biotype"))
 print( Heatmap(asinh(module.grey),column_title=paste0("Grey Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.grey)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.grey)],name="gene_biotype"))
 print(Heatmap(asinh(module.turquoise),column_title=paste0("Turquoise Module Repeat Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=6))+Heatmap(rowRanges(rexp)$tx_biotype[rowRanges(rexp)$tx_id%in%rownames(module.turquoise)],name="tx_biotype")+ Heatmap(rowRanges(rexp)$gene_biotype[rowRanges(rexp)$tx_id%in%rownames(module.turquoise)],name="gene_biotype"))
dev.off()


} ##main
