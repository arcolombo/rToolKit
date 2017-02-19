#' @title this analyzes correlations of repeat and gene modules of fixed interest
#' @description after running wgcna and wrcna we are interested in finding associaitons between one fixed gene with all the repeat modules.  given a fixed module of interest we correlate and test for comparability between every repeat module.  If there is no contrast between the gene modules, i.e. qusage contrasts are not present, then we use the wgcna_goEnrich caller
#' @param kexp this should be the full stage kexp, could even be a tcga kexp
#' @param lnames is the wgcna.R output and rnames is the wrcna.R call
#' @param wgcnaDbName this is the sql db of the wgcna
#' @param wrcnaDbName sql db of the wrcna repeat db
#' @param qusageDbName this is the sql db of a qusage db. not required if there is a missing qusage db one can call goEnrich for preliminary analyses
#' @param geneModules this should be single gene module of interest.  we can examine one gene module's association to all the repeat modules.  so geneModule should be a single color character
#' @param repeat.cutoff the cpm floor threshold if the cut off filters out elements that were not filtered in the wgcna method call, this will fail.  you must match the read.cutoff in the wgcna data base.
#' @param how  either cpm or tpm, tpm recommended for heatmap visualizations
#' @param enrichmentCaller either go or qusage, for some cases go enrichment is an interest to study
#' @param comparison1 this is the first comparison group the comparison key is in lower case letters 'mds' etc.
#' @param comparison2 this is the second comparison group, if there is only 1 comparison group, then leave this as "blast"; the method will query the 'blast' contrast and come up with an empty set and discard it.
#' @param openDevice boolean if true then the images are opened on a device along with printing to pdf.  should be false for exploring all colors
#' @param p.value  an interger significance alpha level
#' @import WGCNA
#' @import edgeR
#' @import pvclust
#' @import ComplexHeatmap
#' @import dendsort
#' @export
singleModuleAnalysis<-function(kexp,lnames,rnames,wgcnaDbName="wgcna.blast.DbLite.sqlite",wrcnaDbName="wrcna.blast.DbLite.cpm.sqlite",qusageDbName="qusage.blast.DbLite.cpm.sqlite", geneModules=NULL,read.cutoff=4,how="tpm",enrichmentCaller=c("go","qusage"),comparison1="phsc",comparison2="blast",openDevice=FALSE,p.value=0.1,design=metadata(kexp)$design,hubPath="~/Documents/Arkas-Paper-Data/AML-bonemarrow-LSCs/patient-plot-data/wgcna_data/repeat_short_list_biotypes/Buen_TCGA_Normal_Karytype_Analysis_ens84_rep2103/GO-analysis/BuenRostro_GO_analysis/Re-Analysis-Blast-LSC/",hubFile="top_65_HubGenes.AML.csv" ){
  ##FIX ME: generically handle repeat modules of various sizes.
  
repeatColors<-listModuleColors(wgcnaDbLite(wrcnaDbName))
repeatModules<-repeatColors
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
  key.id<-which(paste0("ME",geneModules)==colnames(MEs))  
  
  gene.color<-modulesBy(wgcnaDbLite(wgcnaDbName),Module.color=geneModules)
 ##print out the csv of the gene correlation matrix for each module
  repeat.module.list<-lapply(repeatModules,function(x) modulesBy(wgcnaDbLite(wrcnaDbName),Module.color=x))
  names(repeat.module.list)<-repeatModules
  ##print out the csv of the repeat correlation matrix 

 
 ##print out the moduleTriat cor for the specific repeat module this is the key.
   repeat.key<-renameRepeatModuleColors(rnames=rnames,rdbName=wrcnaDbName,wgcnaDbName=wgcnaDbName,geneModules=geneModules,MEs=MEs)


  if(enrichmentCaller=="qusage"){
  ## print enrichment data to screen and write out to CSV
   pHSC<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast=comparison1,p.value=p.value)
  if(nrow(pHSC)>0){
  pHSC<-data.frame(pHSC,contrast=comparison1,color=geneModules)
    }
  blast<-pathways(qusageDbLite(qusageDbName),Module.color=geneModules,tx.Biotype=geneModules,contrast=comparison2,p.value=p.value)  
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
  write.csv(pathWays ,file=paste0(geneModules,"_",key.id,"_",how,"_moduleWiseAnalysisPathways.csv"))
  if(openDevice==TRUE){
   readkey()
   }
  }else{
  #if no qusage contrasts, call goEnrich and print out the fixed color module
  goenriched<-goEnrich(wgcnaDbLite(wgcnaDbName),Module.color=geneModules)
  print(goenriched)
  write.csv(as.data.frame(goenriched) ,file=paste0(geneModules,"_",key.id,"_",how,"_GOmoduleWiseAnalysis.csv"))
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
  ##write out the MEs data used in bar plotting
  write.csv(MEs,file=paste0("All.Gene.ModuleEigenValues.perSample.expression.summarization.csv"))
 rMEs2<-rMEs
  rownames(rMEs2)<-rownames(MEs)
  write.csv(rMEs2,file="All.Repeat.ModuleEigenValues.perSample.repeat.expression.summary.csv")

  color.ID<-which(paste0("ME",geneModules)==colnames(MEs))
 ### barplot of sample module eigen values for a given module
  if(openDevice==TRUE){
  dev.new()
  par(mar=c(4,10,4,4))
  barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.8,las=1,main=paste0(geneModules," (",key.id,") Gene Module Eigenvalues"),space=1)

   
 }
  readkey()

  
  ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
 # if(openDevice==TRUE){
  ###heatmap of module genes
 dev.new()
 ####Heatmap DE genes in the module 


  m.kexp<- kexp[rowRanges(kexp)$gene_id%in%rownames(module.cpm),]

   tpm<-collapseTpm(m.kexp,"gene_name")


  mod.tpm<-tpm
 rpt.pv<-pvclust(log(1+mod.tpm),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+mod.tpm))),isReverse=TRUE)
  rh.rpt<-Heatmap((mod.tpm),
                  name="tpm",
                  cluster_rows=rpt.dend,
                  cluster_columns=rpt.pv$hclust,
                  column_title=paste0("Gene Module ",key.id," (",geneModules,") "," Full Module Expression"),
                  column_title_gp=gpar(fontsize=15),
                 row_names_gp=gpar(fontsize=7),
                 column_names_gp=gpar(fontsize=13),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))
  

  gwa<-geneWiseAnalysis(m.kexp,design=metadata(m.kexp)$design,how="tpm",read.cutoff=2,species="Homo.sapiens",adjustBy="none")
  de.tpm<-tpm[rownames(tpm)%in%as.character(gwa$limmaWithMeta[,8]),]
  top.tpm<-de.tpm
  if(nrow(top.tpm)>15){
   de.pv<-pvclust(log(1+top.tpm),nboot=100)
  de.dend<-dendsort(hclust(dist(log(1+top.tpm))),isReverse=TRUE)
  de.rpt<-Heatmap(log(1+top.tpm),
                  name="tpm",
                  cluster_rows=de.dend,
                  cluster_columns=de.pv$hclust,
                  column_title=paste0("Gene Module ",key.id," (",geneModules,") "," DE Gene Expression"),
                  column_title_gp=gpar(fontsize=15),
                 row_names_gp=gpar(fontsize=8),
                 column_names_gp=gpar(fontsize=13),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))

   }else{
 de.rpt<-Heatmap(log(1+top.tpm),
                  name="tpm",
                  column_title=paste0("Gene Module ",key.id," (",geneModules,") "," DE Gene Expression"),
                  column_title_gp=gpar(fontsize=15),
                 row_names_gp=gpar(fontsize=8),
                 column_names_gp=gpar(fontsize=13),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))


}


 
 

  draw(rh.rpt)
  readkey()
  draw(de.rpt)
   
  ##grab hubs

hub<-read.csv(paste0(hubPath,"/",hubFile),header=FALSE)
pale.id<-which(geneModules==as.character(hub$V1))
pale<-as.character(hub[pale.id:(pale.id+66),4])
pale<-pale[!is.na(pale)]

  hub.tpm<-tpm[rownames(tpm)%in%pale,]

  hub.pv<-pvclust(log(1+hub.tpm),nboot=100)
  hub.dend<-dendsort(hclust(dist(log(1+hub.tpm))),isReverse=TRUE)
  hub.rpt<-Heatmap(log(1+hub.tpm),
                  name="tpm",
                  cluster_rows=hub.dend,
                  cluster_columns=hub.pv$hclust,
                  column_title=paste0("Module Hub ",key.id," (",geneModules,") "," Hub Gene Expression"),
                  column_title_gp=gpar(fontsize=15),
                 row_names_gp=gpar(fontsize=8),
                 column_names_gp=gpar(fontsize=13),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))

readkey()
 draw(hub.rpt)



 
 ### PRINT TO PDF########3
  pdf(paste0(geneModules,"_",key.id,"_",how,"_",enrichmentCaller,"_EigenValueAnalysis.pdf"))
  par(mar=c(4,10,4,4))
    barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.8,las=1,main=paste0(geneModules," (",key.id,") Gene Module Eigenvalues"),space=1)
  dev.off()
 
  pdf(paste0(geneModules,"_",key.id,"_",how,"_",enrichmentCaller,"_ModuleExpressions.pdf"),height=11)

 print(rh.rpt)
 print(de.rpt)
 print(hub.rpt)
  dev.off() 
##heatmap of repeat modules
   

} ##main


