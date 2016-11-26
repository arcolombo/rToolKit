#' @title produces a correlation map given wgcna network data using a different library ComplexHeatmap
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image with the columns ordered.
#' @param lnames this is the results from wgcna.R
#' @param read.cutoff  integer for min cutoff
#' @param plotDot boolean this plots the median correlation score for each biotypes in a line plot
#' @param recalc boolean if truen then will recalculate the bicor and pvalues
#' @param targets this is a character vector of modules of interest to subset the corMap
#' @param orderBiotype either Alu or L1.
#' @import WGCNA
#' @import ComplexHeatmap
#' @import pvclust
#' @import dendsort
#' @export
#' @return images of eigengenes
wgcna_Heatcor<-function(lnames=NULL,read.cutoff=2,recalc=FALSE,targets=NULL,how=how,pathwaysToPick=c("immune","inflam","apopto","death","kappab","wound"),pathPairing=c(1,1,2,2,3,4),qdbname=NULL){
##FIX ME:  add a star for pvalues less than 0.05 in the cell_function Heat
stopifnot(length(pathwaysToPick)==length(pathPairing))
 
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
  if(is.null(targets)==FALSE) {
  moduleTraitCor<-moduleTraitCor[rownames(moduleTraitCor)%in%targets,]
  moduleTraitPvalue<-moduleTraitPvalue[rownames(moduleTraitPvalue)%in%targets,]
 }
 # if(orderBiotype=="L1"){
 #moduleTraitCor<-moduleTraitCor[order(moduleTraitCor[,grep("L1",colnames(moduleTraitCor))],decreasing=TRUE),]
# } else if(orderBiotype=="Alu"){
#  moduleTraitCor<-moduleTraitCor[order(moduleTraitCor[,grep("Alu",colnames(moduleTraitCor))],decreasing=TRUE),]
# }

 moduleTraitPvalue.id<-match(rownames(moduleTraitCor),rownames(moduleTraitPvalue))
 moduleTraitPvalue<-moduleTraitPvalue[moduleTraitPvalue.id,]
 
 ###label the rows of the heatmap by numbers
key<-data.frame(module=rownames(moduleTraitCor),id=seq(1:length(rownames(moduleTraitCor))),stringsAsFactors=FALSE)

 ###color code the rownames based on pathPairing
 ##add a legened
 ##create  a sub plot that shows the pathwaysToPick specific module Heatmap Correlations with altered significance level.
xN<-list()
for(i in 1:length(pathwaysToPick)){
  xnam<-names(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i]))
  if(length(xnam)>0){
  xnam<-xnam[sapply(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i]),function(x) nrow(x)>=1)]
  
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
 modA<-HeatmapAnnotation(df=df_annot,which="row")
####
  moduleTraitCor2<-moduleTraitCor
  key.id2<-match(rownames(moduleTraitCor2),key$module)
  rownames(moduleTraitCor2)<-key$id[key.id2]
  par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
  plot.new()
 if(nrow(moduleTraitCor2)>6){
  x.pv<-pvclust(moduleTraitCor2,nboot=100)
  
  heatCor<-Heatmap(moduleTraitCor2,cluster_columns=x.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
  }
  readkey()

#####add Heatcor of row annotation subset. 

toGo<-data.frame(module=unique(unlist(xN)),stringsAsFactors=FALSE)
 mC<-moduleTraitCor[rownames(moduleTraitCor)%in% paste0("ME",toGo$module),]
 mc.id<-match(rownames(mC),key$module)
 rownames(mC)<-key$id[mc.id]


df_annot2<-data.frame(module=rownames(moduleTraitCor))
rownames(df_annot2)<-rownames(moduleTraitCor)
for(i in 1:length(names(xN))){

df_annot2<-cbind(df_annot2,-1)
colnames(df_annot2)[i+1]<-names(xN)[i]

 df_annot2[ rownames(df_annot2)%in% paste0("ME", xN[[i]]  )  ,i+1  ]<-1
}
 df_annot2$module<-NULL
 key.id1<-match(rownames(df_annot2),key$module)
 rownames(df_annot2)<-key$id[key.id1]
 df_annot2<-df_annot2[rownames(df_annot2)%in%rownames(mC),]
 modB<-HeatmapAnnotation(df=df_annot2,which="row")

###add activation direction
activation.direction<-matrix(data=0,nrow=nrow(mC),ncol=length(names(xN)))
rownames(activation.direction)<-rownames(mC)
colnames(activation.direction)<-names(xN)
for(i in 1:length(pathwaysToPick)) {
 df<-pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i])
  if(length(df)==0){
  next
 }
 df<-df[sapply(df,function(x) nrow(x)>=1)]

 activation.df<-sapply(df,function(x) as.numeric(x[grep("Total",rownames(x)),]$logFC)/as.numeric(x[grep("Total",rownames(x)),]$query.size))
  names(activation.df)<-paste0("ME",names(activation.df))
 dir.id<-match(names(activation.df),key$module)
  names(activation.df)<-dir.id
 dir.id2<-match(names(activation.df),rownames(activation.direction))

 activation.direction[dir.id2, which(colnames(activation.direction)==pathwaysToPick[i])]<-activation.df
}

########
  activation.Heat<-Heatmap(asinh(activation.direction),col=colorRamp2(c(-2,0,2),c("black","white","orange")),name="Pathway logFC",column_title="Pathway Activity"  )
  par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
  plot.new()
  if(nrow(mC)>4){
  x.pv2<-pvclust(mC,nboot=200)

  heatCor2<-Heatmap(mC,cluster_columns=x.pv2$hclust,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Immune-Related ",how," Biotype relationships"))
  print(modB+heatCor2+activation.Heat)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Immune-Related ",how," Biotype relationships"))
  print(modB+heatCor2+activation.Heat)
  }
  readkey()
#######

### adjust p.values 
### change the rownames of the heatmap to 1:n numerics as the last step. 


 
  pdf(paste0("Heat_correlation_",how,"_plots.pdf"),width=12,height=12)
    par(mar = c(6, 10, 3, 3));
     print(modA+heatCor)
    print(modB+heatCor2+activation.Heat)
   dev.off()
 



} ###main



