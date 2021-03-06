#' @title produces a correlation map given wgcna network data using a different library ComplexHeatmap
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image with the columns ordered. This will produce the HeatCorrelation association table with the pathway information alongside.  this will also produce the repeat module renaming call to print the summary of repeat module tables as a repeat correlation summary.
#' @param lnames this is the results from wgcna.R
#' @param read.cutoff  integer for min cutoff
#' @param recalc boolean if truen then will recalculate the bicor and pvalues
#' @param how tpm or cpm, this is used to print the heatmap, tpm is better
#' @param pathwaysToPick this will query the qusageDbLite database for the key pathway character name signature
#' @param pathPairing vector if pathways are related then they should have matching pairing integers for instance cell death and apoptosis pathways could be paired.
#' @param dbname  the gene module wgcnaDbLite sqlite database
#' @param qdbname  the qusageDbLite sqlite database name
#' @param rdbname the repeat module wrcnaDbLite sqlite database
#' @param averagePathwayFC is a boolean, where you wish to average pathways related to a class. this should be true if the canonical merged pathways is not used. the canonical merged pathways should not use averaging because they combine all related pathways.  if the c5.all was used in the qdb creation, set averaging to TRUE.
#' @import WGCNA
#' @import ComplexHeatmap
#' @import pvclust
#' @import dendsort
#' @import plyr
#' @import circlize
#' @export
#' @return images of eigengenes
wgcna_HeatCorFinal<-function(lnames=NULL,rnames=NULL, read.cutoff=1,recalc=FALSE,how=how,pathwaysToPick=c("immune_response","inflam","apopto","death","kappab","wound"),pathPairing=c(1,1,2,2,3,4),dbname=NULL,qdbname=NULL,rdbname=NULL,p.value=1,minRank=0,averagePathwayFC=FALSE){
##FIX ME:  add a star for pvalues less than 0.05 in the cell_function Heat
stopifnot(length(pathwaysToPick)==length(pathPairing))
 # require(plyr)
 # require(ComplexHeatmap)
 # require(circlize)
 # require(pvclust) 
 # require(WGCNA)
if(is.null(lnames)==TRUE){
cat("Please call wgcna.R prior calling a heatmap...\n")
}
###declare needed objects from load
if(is.null(lnames)==FALSE){
message(paste0("found lnames"))
} else {
 stop("did not find lnames, please run wgnca")
}
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
# Will display correlations and their p-values
# Define numbers of genes and samples
  nGenes= ncol(datExpr);
  nGenes= ncol(datExpr);
  nSamples = nrow(datExpr);
 # Recalculate MEs with color labels
  if(recalc==TRUE){
   moduleTraitCor = bicor(MEs, datTraits, use = "all.obs");
   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
   } else{
   moduleTraitCor<-lnames[["moduleTraitCor"]]
   moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
   }

 moduleTraitPvalue.id<-match(rownames(moduleTraitCor),rownames(moduleTraitPvalue))
 moduleTraitPvalue<-moduleTraitPvalue[moduleTraitPvalue.id,]
 write.csv(moduleTraitPvalue,file="Correlation_Matrix_PValues_ModuleRepeat-CPM-Biotype-TxBiotype.csv")

 ###label the rows of the heatmap by numbers
    key<-data.frame(module=rownames(moduleTraitCor),id=seq(1:length(rownames(moduleTraitCor))),stringsAsFactors=FALSE)
 write.csv(key,file="ModuleNameColor_key_ModuleNumber.csv")

 ###color code the rownames based on pathPairing
 ##add a legened
 ##create  a sub plot that shows the pathwaysToPick specific module Heatmap Correlations with altered significance level.
  xN<-list()
  for(i in 1:length(pathwaysToPick)){
    xnam<-names(pickPathway(qusageDbLite(qdbname),
                keyWord=pathwaysToPick[i],
                p.value=p.value))
  if(length(xnam)>0){
  xnam<-xnam[sapply(pickPathway(qusageDbLite(qdbname),
                    keyWord=pathwaysToPick[i],
                    p.value=p.value),function(x) median(x[which(rownames(x) !="Total"),"ranking"])>minRank)    ]
     xN[[i]]<-xnam
     names(xN)[i]<-pathwaysToPick[i]
     xN<-xN[sapply(xN,function(x) length(x)>0)]
  }
}
###row Annotation left side
  df_annot<-data.frame(module=rownames(moduleTraitCor))
  rownames(df_annot)<-rownames(moduleTraitCor)
  for(i in 1:length(names(xN))){
  df_annot<-cbind(df_annot,-1)
  colnames(df_annot)[i+1]<-names(xN)[i]
  df_annot[ rownames(df_annot)%in% paste0("ME", xN[[i]]  )  ,i+1  ]<-1
 }
   df_annot$module<-NULL
   key.id1<-match(rownames(df_annot),key$module)
   rownames(df_annot)<-key$id[key.id1]
   modA<-HeatmapAnnotation(df=df_annot,which="row" )
####
  
  moduleTraitCor2<-moduleTraitCor
  key.id2<-match(rownames(moduleTraitCor2),key$module)
  rownames(moduleTraitCor2)<-key$id[key.id2]
 write.csv(moduleTraitCor2,file="Correlation_Matrix_ModuleRepeat-CPM-Biotype-TxBiotype.csv")

  par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
   plot.new()
   if(nrow(moduleTraitCor2)>6){
   x.pv<-pvclust(moduleTraitCor2,nboot=100)
   colnames(moduleTraitCor2)<-gsub("Repetitive element","Rptv. Element",colnames(moduleTraitCor2))
  colnames(moduleTraitCor2)<-gsub("Endogenous Retrovirus","Endg. Retrovirus",colnames(moduleTraitCor2))
  heatCor<-Heatmap(moduleTraitCor2,column_names_gp=gpar(fontsize=10),cluster_columns=x.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="correlation(x)", 
   heatmap_legend_param=list(color_bar="continuous"),
   column_title = paste0("Module-Repeat ",how," Biotype relationships (*<0.05)"),
        cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(moduleTraitCor2,nSamples)
                   if(weighted.Pvalue[i,j]<0.05){
                    grid.text("*",x,y-unit(0.005,'npc') )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }     
         )
  print(modA+heatCor)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
  }
  readkey()
###add activation direction
 mC<-moduleTraitCor
  activation.direction<-matrix(data=0,nrow=nrow(moduleTraitCor),ncol=length(names(xN)))
  rownames(activation.direction)<-rownames(mC)
  colnames(activation.direction)<-names(xN)
  rownames(activation.direction)<-key$id[match(key$module,rownames(activation.direction))]
  for(i in 1:length(pathwaysToPick)) {
  df<-pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i],p.value=p.value)
  df<-df[which(names(df)!="kexp")]
  ##print to csv
  activation.output<-ldply(df,data.frame)
  write.csv(activation.output,file=paste0("Activation.Direction.Pathways.",pathwaysToPick[i],"_AcrossModules.csv"  ) )
  if(length(df)==0){
  next
 }
 df<-df[sapply(df,function(x) median(x[which(rownames(x) !="Total"),"ranking"])>minRank)  ]
  ###average pathways is c5.all.symbols.gmt was used.  if the canonical merged GO.gmt was used do not average.
 if(averagePathwayFC==TRUE){
    df.path<-lapply(df,function(x) x[which(rownames(x) !="Total"),] )
    sig.df<-lapply(df.path,function(x) x[which(x[,"pvalue"]<=p.value),])
   activation.df<-sapply(sig.df,function(x) sum(x$logFC))
  }else{
  merged.df<-lapply(df,function(x) x["Total",])
 sig.df<-merged.df[sapply(merged.df,function(x) x[,"pvalue"]<=0.05)]
 activation.df<-sapply(sig.df,function(x) as.numeric(x[grep("Total",rownames(x)),]$logFC))
   }
 names(activation.df)<-paste0("ME",names(activation.df))
 dir.id<-match(names(activation.df),key$module)
  names(activation.df)<-dir.id
 dir.id2<-match(names(activation.df),rownames(activation.direction))

 activation.direction[dir.id2, which(colnames(activation.direction)==pathwaysToPick[i])]<-activation.df
}

########
  activation.Heat<-Heatmap(asinh(activation.direction),col=colorRamp2(c(-2,0,2),c("black","white","orange")),name="Pathway logFC",column_title="Pathway Activity",show_row_names=FALSE  )

    print(heatCor+activation.Heat)
    readkey()
#######

 
  pdf(paste0("Heat_correlation_",how,"_plots.mergedPathwaySignificant.pdf"),width=12,height=12)
    par(mar = c(9, 10, 3, 3));
    print(heatCor+activation.Heat)
   dev.off()
 



} ###main



