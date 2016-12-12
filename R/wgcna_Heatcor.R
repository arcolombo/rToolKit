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
#' @import WGCNA
#' @import ComplexHeatmap
#' @import pvclust
#' @import dendsort
#' @import plyr
#' @import circlize
#' @export
#' @return images of eigengenes
wgcna_Heatcor<-function(lnames=NULL,rnames=NULL, read.cutoff=1,recalc=FALSE,how=how,pathwaysToPick=c("immune","inflam","apopto","death","kappab","wound"),pathPairing=c(1,1,2,2,3,4),dbname=NULL,qdbname=NULL,rdbname=NULL){
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
  xnam<-names(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i]))
  if(length(xnam)>0){
  xnam<-xnam[sapply(pickPathway(qusageDbLite(qdbname),keyWord=pathwaysToPick[i]),function(x) median(x[which(rownames(x) !="Total"),"ranking"])>40)    ]
  
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
  heatCor<-Heatmap(moduleTraitCor2,column_names_gp=gpar(fontsize=10),cluster_columns=x.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships (*<0.06)"),cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(moduleTraitCor2,nSamples)
                   if(weighted.Pvalue[i,j]<0.06){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
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

#####add Heatcor of row annotation subset. 

  toGo<-data.frame(module=unique(unlist(xN)),stringsAsFactors=FALSE)
  mC<-moduleTraitCor[rownames(moduleTraitCor)%in% paste0("ME",toGo$module),]
  mc.id<-match(rownames(mC),key$module)
  rownames(mC)<-key$id[mc.id]
 colnames(mC)<-gsub("Repetitive element","Rptv. Element",colnames(mC))
  colnames(mC)<-gsub("Endogenous Retrovirus","Endg. Retrovirus",colnames(mC))
 write.csv(mC,file="Module-Repeat-CPM-Biotype.relationships.csv")
 weighted.Pvalue<-corPvalueStudent(mC,nSamples)
  write.csv(weighted.Pvalue,"Module-Repeat-CPM-Biotype-releaitonships.P.value.csv")
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
  df<-df[which(names(df)!="kexp")]
  ##print to csv
  activation.output<-ldply(df,data.frame)
  write.csv(activation.output,file=paste0("Activation.Direction.Pathways.",pathwaysToPick[i],"_AcrossModules.csv"  ) )
  if(length(df)==0){
  next
 }
 df<-df[sapply(df,function(x) median(x[which(rownames(x) !="Total"),"ranking"])>40)  ]

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

  heatCor2<-Heatmap(mC,cluster_columns=x.pv2$hclust,cluster_rows=FALSE,column_names_gp=gpar(fontsize=10),row_names_side="left",name="correlation(x)", column_title = paste0("Immune-Related ",how," (*<0.06)"), cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(mC,nSamples)
                   if(weighted.Pvalue[i,j]<0.06){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }

       )
  print(modB+heatCor2+activation.Heat)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Immune-Related ",how," Biotype relationships"))
  print(modB+heatCor2+activation.Heat)
  }
  readkey()
#######

 ##traitCorRenamed is the slot for custom renaming of module repeat identity of family
 if(is.null(rnames[["traitCorRenamed"]])==FALSE){
  rTraitCor<-rnames[["traitCorRenamed"]]
  }else{
  ##FIX ME: if traitCorRenamed is not found call renameScript
  repeat.key<-renameRepeatModuleColors(rnames=rnames,rdbName=rdbname,wgcnaDbName=dbname,geneModules="black",MEs=MEs)
  rTraitCor<-rnames[["moduleTraitCor"]]
  rename.id<-match(rownames(rTraitCor),rownames(repeat.key))
   rownames(rTraitCor)<-repeat.key$leading.terms[rename.id]
   write.csv(rTraitCor,file="repeat.biotype.correlations.renamedTo.leading.biotypeCorrelations.csv")
  }
   if(ncol(as.data.frame(moduleTraitCor))>1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%colnames(moduleTraitCor)]
  }else if(ncol(as.data.frame(moduleTraitCor))==1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%names(moduleTraitCor)]
  }
  write.csv(rTraitCor,file="Correlation-RepeatModules-RepeatBiotypes-Repeat_Module_Correlations.csv")
  corrMap<-bicor(t(moduleTraitCor),t(rTraitCor))
  corrMap.pvalue<-corPvalueStudent(corrMap,nSamples)
  corrMap.pvalue[which(is.na(corrMap.pvalue))]<-1
  write.csv(corrMap,file="Correlation-GeneModules-RepeatModules.csv")
  write.csv(corrMap.pvalue,file="Correlation-GeneModules-RepeatModules.P.values.csv")

 ####
  corrMap2<-corrMap
  key.id2<-match(rownames(corrMap2),key$module)
  rownames(corrMap2)<-key$id[key.id2]
    write.csv(corrMap2,file="Correlation-GeneModules-RepeatModules-Pvalues.csv")

  par(mar = c(10, 10, 3, 3));
 # Display the correlation values within a heatmap plot
  plot.new()
  cor.pv<-pvclust(corrMap2,nboot=100)
 map.heatCor<-Heatmap(corrMap2,cluster_columns=cor.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="correlation(x)",column_names_gp=gpar(fontsize=7.5), column_title = paste0("Repeat Level ",how," Module (*<0.06)"),cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(corrMap2,nSamples)
                   if(weighted.Pvalue[i,j]<0.06){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }


  )
  print(modA+map.heatCor)

###########

##summary module repeat 
 corrMap.sub<-corrMap2
 key.id3<-match(rownames(mC),rownames(corrMap.sub))
 corrMap.sub<-corrMap.sub[key.id3,]
 if(nrow(corrMap.sub)>4){
  sub.pv<-pvclust(corrMap.sub,nboot=200)

  sub.heatCor2<-Heatmap(corrMap.sub,cluster_columns=sub.pv$hclust,cluster_rows=FALSE,row_names_side="left",name="correlation(x)",column_names_gp=gpar(fontsize=7), column_title = paste0("Immune-Related ",how," Module relationships"), 
cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(corrMap.sub,nSamples)
                   if(weighted.Pvalue[i,j]<0.06){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }

   )
  print(modB+sub.heatCor2+activation.Heat)
  } else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Immune-Related ",how," Biotype relationships"),
 cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(moduleTraitCor,nSamples)
                   if(weighted.Pvalue[i,j]<0.06){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }

 )
  print(modB+sub.heatCor2+activation.Heat)
  }
  readkey()
#############


### adjust p.values 
### change the rownames of the heatmap to 1:n numerics as the last step. 


 
  pdf(paste0("Heat_correlation_",how,"_plots.pdf"),width=12,height=12)
    par(mar = c(9, 10, 3, 3));
     print(modA+heatCor)
    print(modB+heatCor2+activation.Heat)
    print(modA+map.heatCor)
    print(modB+sub.heatCor2+activation.Heat)
   dev.off()
 



} ###main



