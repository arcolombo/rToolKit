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
wgcna_Heatcor<-function(lnames=NULL,read.cutoff=2,plotDot=FALSE,recalc=FALSE,targets=NULL,how=how,reOrder=TRUE,pathwaysToPick=c("immune","inflammation","apoptotic","death","kappab","wound"),pathPairing=c(1,1,2,2,3,4),qdbname=NULL){

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


 ###color code the rownames based on pathPairing
 ##add a legened
 ##create  a sub plot that shows the pathwaysToPick specific module Heatmap Correlations with altered significance level.
xN<-list()
for(i in 1:length(pathwaysToPick)){
xnam<-names(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i])[sapply(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i]),function(x) nrow(x)>2)])
xN[[i]]<-xnam
names(xN)[i]<-pathwaysToPick[i]
}
xN<-xN[sapply(xN,function(x) length(x)>0)]
###row Annotation left side
df_annot<-data.frame(module=rownames(moduleTraitCor))
rownames(df_annot)<-rownames(moduleTraitCor)
for(i in 1:length(names(xN))){

df_annot<-cbind(df_annot,-1)
colnames(df_annot)[i+1]<-names(xN)[i]

 df_annot[ rownames(df_annot)%in% paste0("ME", xN[[i]]  )  ,i+1  ]<-1
}
 df_annot$module<-NULL
 modA<-HeatmapAnnotation(df=df_annot,which="row")
####


  par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
  plot.new()
  if(reOrder==TRUE){
  x.pv<-pvclust(moduleTraitCor,nboot=100)
  
  heatCor<-Heatmap(moduleTraitCor,cluster_columns=x.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
  }
  readkey()

#####add Heatcor of row annotation subset.  
### adjust p.values 
### change the rownames of the heatmap to 1:n numerics as the last step. 


  colorDF<-sapply(MEs,function(x) median(x))
  colorDF<-colorDF[order(colorDF,decreasing=TRUE)]
  if(plotDot==TRUE){
  plot.new()
  par(mar = c(11.5, 6, 3, 3));
  plot(colorDF,main="Median Correlation Per Module")
  axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
  readkey()
  } else {
  pdf(paste0("Heat_correlation_",how,"_plots.pdf"),width=12,height=12)
    par(mar = c(6, 10, 3, 3));
    if(reOrder==TRUE){
     x.pv<-pvclust(moduleTraitCor,nboot=100)
     heatCor<-Heatmap(moduleTraitCor,cluster_columns=x.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(modA+heatCor)
   } else{
   heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))   
    print(modA+heatCor)
  }
  plot(colorDF,main="Median Correlation Per Module")
  axis(1,at=1:length(colorDF),labels=names(colorDF),las=2)
  dev.off()
 }



} ###main



