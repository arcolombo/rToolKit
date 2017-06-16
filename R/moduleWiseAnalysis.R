#' @title this analyzes correlations of repeat and gene modules of fixed interest
#' @description after running wgcna and wrcna we are interested in finding associaitons between one fixed gene with all the repeat modules.  given a fixed module of interest we correlate and test for comparability between every repeat module.  If there is no contrast between the gene modules, i.e. qusage contrasts are not present, then we use the wgcna_goEnrich caller
#' @param kexp this should be the full stage kexp, could even be a tcga kexp
#' @param lnames is the wgcna.R output and rnames is the wrcna.R call
#' @param wgcnaDbName this is the sql db of the wgcna
#' @param wrcnaDbName sql db of the wrcna repeat db
#' @param qusageDbName this is the sql db of a qusage db. not required if there is a missing qusage db one can call goEnrich for preliminary analyses
#' @param geneModules this should be single gene module of interest.  we can examine one gene module's association to all the repeat modules.  so geneModule should be a single color character
#' @param repeatModules this should be the character of all repeat modules
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
moduleWiseAnalysis<-function(kexp,lnames,rnames,wgcnaDbName="wgcnaDbLite.cpm.sqlite",wrcnaDbName="wrcnaDbLite.cpm.sqlite",qusageDbName="qusageDbLite.cpm.sqlite", geneModules=c("darkturquoise","red","grey60","pink","darkred","royalblue","green","magenta","brown"),repeatModules=c("blue","brown","grey","turquoise"),read.cutoff=4,how="tpm",enrichmentCaller=c("go","qusage"),comparison1="phsc",comparison2="blast",openDevice=FALSE,p.value=0.1,design=metadata(kexp)$design,includeRepeatHeatmaps=FALSE ){
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
  write.csv(pathWays ,file=paste0(geneModules,"_",key.id,"_",how,"_moduleWiseAnalysis.csv"))
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

  for(colR in repeatModules){
  dev.new()
    par(mar=c(4,10,4,4),cex.main=0.9)
   color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  ME.id<-which(colnames(rMEs)==(paste0("ME",colR)))
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
  barplot(rMEs[,ME.id],horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0(repeat.key[which(rownames(repeat.key)==paste0("ME",colR)),]," Module correlation ~(",key.id,") ",signif(bicor(rMEs[,ME.id],MEs[,color.ID]),3 )) ,space=1)
   
 }
  readkey()

  }
  ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
 # if(openDevice==TRUE){
  ###heatmap of module genes
 dev.new()
 ####Heatmap DE genes in the module 


  m.kexp<- kexp[rowRanges(kexp)$gene_id%in%rownames(module.cpm),]
  res<-fitBundles(m.kexp,design,bundleID="gene_name",read.cutoff=read.cutoff)
   tpm<-collapseTpm(m.kexp,"gene_name")
   
    res$top <- with(res, topTable(fit, coef=2, p=0.05,adjust.method="none", n=nrow(kexp)))
 if(nrow(res$top)>0 && nrow(res$top)<=64 ){
    top.tpm<-tpm[rownames(tpm)%in%rownames(res$top),]
  }else if(nrow(res$top)>65){
   res$top <- with(res, topTable(fit, coef=2, p=0.05,adjust.method="BH", n=100))
   top.tpm<-tpm[rownames(tpm)%in%rownames(res$top),]
   }else{
  res$top<-byMad(tpm,nrow(tpm)/(1.5))
  top.tpm<-tpm[rownames(tpm)%in%rownames(res$top),]
   }
 # }

  ########
  if(nrow(res$top)>=10){ 
   rpt.pv<-pvclust(log(1+top.tpm),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+top.tpm))),isReverse=TRUE)
  rh.rpt<-Heatmap(log(1+top.tpm),
                  name="log(1+tpm)",
                  cluster_rows=rpt.dend,
                  cluster_columns=rpt.pv$hclust,
                  column_title=paste0("Gene Module ",key.id," (",geneModules,") "," DE Gene Expression"),
                  column_title_gp=gpar(fontsize=15),
                 row_names_gp=gpar(fontsize=8),
                 column_names_gp=gpar(fontsize=13),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))

  }else{
 rh.rpt<-Heatmap(log(1+top.tpm),
                  name="log(1+tpm)",
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
 # print(Heatmap(asinh(module.cpm),column_title=paste0(geneModules," (",key.id,") Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=5) ) )
 ##heatmap of repeat modules
 # dev.new()

   tx_biotype_color_master<-repeatBiotypeColorPalette(kexp)
    tx.color.master<-tx_biotype_color_master
   gene_biotype_color_master<-repeatFamilyColorPalette(kexp)
   gn.color.master<-gene_biotype_color_master 
  ##call DE repeat expression using drawHeatmap
   for(colR in repeatModules){

  color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,] #all the repeats in the module color
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
   r.title<-repeat.key[grep(paste0("ME",colR),rownames(repeat.key)),]
   mr.kexp<- kexp[rowRanges(kexp)$tx_id%in%rownames(module.color),]
  rps<-fitBundles(mr.kexp,design,bundleID="gene_id",read.cutoff=read.cutoff)
  rps$top <- with(rps, topTable(fit, coef=2, p=0.05,adjust.method="none", n=nrow(kexp)))

   if(nrow(rps$top)>=10){
   print("printing DE repeats")
  colR.map<-drawHeatmap(mr.kexp,tags=rps$top,byType="tpm",title1=paste0(r.title," DE Repeat Module Expression"),tx_biotype_color_master=tx_biotype_color_master,gene_biotype_color_master=gene_biotype_color_master)  
  }else{
  print(paste0(nrow(rps$top),": not enough DE repeats"))
  if(openDevice==TRUE){
  dev.new()
  }
  rpt.pv<-pvclust(log(1+module.color),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+module.color))),isReverse=TRUE)
  non.de.map<-Heatmap(log(1+module.color),
    column_title=paste0(repeat.key[which(rownames(repeat.key)==paste0("ME",colR)),]," Module Repeat Expression"),
         name="log(1+x)",
         cluster_rows=rpt.dend,
         cluster_columns=rpt.pv$hclust,
         row_names_gp=gpar(fontsize=8),
          column_names_gp=gpar(fontsize=10),
           column_title_gp=gpar(fontsize=9))

 txb.df<-data.frame(transcript_biotype=rowRanges(kexp)[rownames(module.color)]$tx_biotype)
 txb.sel.col=list(transcript_biotype=tx.color.master)
 id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(module.color)]$tx_biotype
 rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))

 gn.df<-data.frame(gene_biotype=rowRanges(kexp)[rownames(module.color)]$gene_biotype)
 gn.sel.col=list(gene_biotype=gn.color.master)
 id3<-names(gn.sel.col[[1]])%in%rowRanges(kexp)[rownames(module.color)]$gene_biotype
 rA2<-rowAnnotation(gn.df,col=list(gene_biotype=(gn.sel.col[[1]][id3])))
 if(openDevice==TRUE){ 
 draw(non.de.map+rA+rA2)
 readkey()
  }
      }
    }##colR loop
# }

 
 ### PRINT TO PDF########3
  pdf(paste0(geneModules,"_",key.id,"_",how,"_",enrichmentCaller,"_EigenValueAnalysis.pdf"))
  par(mar=c(4,10,4,4))
    barplot(MEs[,color.ID],horiz=T,names.arg=rownames(MEs),cex.names=0.8,las=1,main=paste0(geneModules," (",key.id,") Gene Module Eigenvalues"),space=1)

    for(colR in repeatModules){
    par(mar=c(4,10,4,4),cex.main=0.9)
   color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,]
  ME.id<-which(colnames(rMEs)==(paste0("ME",colR)))
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
   barplot(rMEs[,ME.id],horiz=T,names.arg=colnames(kexp),cex.names=0.8,las=1,main=paste0(repeat.key[which(rownames(repeat.key)==paste0("ME",colR)),]," Module correlation ~(",key.id,") ",signif(bicor(rMEs[,ME.id],MEs[,color.ID]),3 )) ,space=1)
   }
  dev.off()
 ##the bar plots of the eigenvalues for repeats that are highly correlated to the gene moduleEigenvalue should have co-expression with similiar covariance 
 # print(Heatmap(asinh(module.cpm),column_title=paste0(geneModules," (",key.id,") Module Gene Expression"),name="asinh(x)",row_names_gp=gpar(fontsize=5) ) )
 
  pdf(paste0(geneModules,"_",key.id,"_",how,"_",enrichmentCaller,"_ModuleExpressions.pdf"),height=11)

 print(rh.rpt)
 ##heatmap of repeat modules
   
  for(colR in repeatModules){
color.id<-match(repeat.module.list[[colR]]$gene_id,rownames(rpm))
  module.color<-rpm[color.id,] #all the repeats in the module color
  stopifnot(all(rownames(module.color)==repeat.module.list[[colR]]$gene_id)==TRUE)
   r.title<-repeat.key[grep(paste0("ME",colR),rownames(repeat.key)),]
   mr.kexp<- kexp[rowRanges(kexp)$tx_id%in%rownames(module.color),]
  rps<-fitBundles(mr.kexp,design,bundleID="gene_id",read.cutoff=read.cutoff)
  rps$top <- with(rps, topTable(fit, coef=2, p=0.05,adjust.method="none", n=nrow(kexp)))

   if(nrow(rps$top)>=10){
  heat.map<-drawHeatmap(mr.kexp,tags=rps$top,byType="tpm",title1=paste0(r.title," DE Repeat  Expression"),tx_biotype_color_master=tx_biotype_color_master,gene_biotype_color_master=gene_biotype_color_master)
  print(heat.map)
  }else{
 # dev.new()
  rpt.pv<-pvclust(log(1+module.color),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+module.color))),isReverse=TRUE)
  non.de.map<-Heatmap(log(1+module.color),
    column_title=paste0(repeat.key[which(rownames(repeat.key)==paste0("ME",colR)),]," Repeat Expression"),
         name="log(1+x)",
         cluster_rows=rpt.dend,
         cluster_columns=rpt.pv$hclust,
         row_names_gp=gpar(fontsize=8),
          column_names_gp=gpar(fontsize=10),
            column_title_gp=gpar(fontsize=9))

 txb.df<-data.frame(transcript_biotype=rowRanges(kexp)[rownames(module.color)]$tx_biotype)
 txb.sel.col=list(transcript_biotype=tx.color.master)
 id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(module.color)]$tx_biotype
 rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))

 gn.df<-data.frame(gene_biotype=rowRanges(kexp)[rownames(module.color)]$gene_biotype)
 gn.sel.col=list(gene_biotype=gn.color.master)
 id3<-names(gn.sel.col[[1]])%in%rowRanges(kexp)[rownames(module.color)]$gene_biotype
 rA2<-rowAnnotation(gn.df,col=list(gene_biotype=(gn.sel.col[[1]][id3])))

 print(non.de.map+rA+rA2)
# readkey()
      }
    }##colR loop
 dev.off()

} ##main


