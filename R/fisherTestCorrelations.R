#' @title tests the correlation and enrichment association from the correlation map 
#' @description after calling CorrelationHeatmap we can test the correlations association to the enrichment activation. This uses a fisher test to test each TE type individually association to the TE type's enrichment activation directions.  finally we aggregate all the TE types count data and perform a global chi-sq test for an overall association.
#' @param lnames this is the results from wgcna.R
#' @param read.cutoff  integer for min cutoff
#' @param recalc boolean if truen then will recalculate the bicor and pvalues
#' @param how tpm or cpm, this is used to print the heatmap, tpm is better
#' @param pathwaysToPick this will query the qusageDbLite database for the key pathway character name signature
#' @param pathPairing vector if pathways are related then they should have matching pairing integers for instance cell death and apoptosis pathways could be paired.
#' @param dbname  the gene module wgcnaDbLite sqlite database
#' @param qdbname  the qusageDbLite sqlite database name
#' @param rdbname the repeat module wrcnaDbLite sqlite database
#' @param p.value this is the p.value threshold for selecting the inferred function of a given module.
#' @param minRank In a given module it can filter based on the competitive rank of where the pathwaysToPick are ranked in a given module.
#' @import WGCNA
#' @import ComplexHeatmap
#' @import pvclust
#' @import dendsort
#' @import plyr
#' @import circlize
#' @export
#' @return images of eigengenes
fisherTestCorrelations<-function(lnames=NULL,rnames=NULL, read.cutoff=1,recalc=FALSE,how=how,pathwaysToPick=c("Immune","Inflammation"),pathPairing=c(1,1),dbname=NULL,qdbname=NULL,rdbname=NULL,p.value=0.1,minRank=0){
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
 #write.csv(moduleTraitPvalue,file="Correlation_Matrix_PValues_ModuleRepeat-CPM-Biotype-TxBiotype.csv")

 ###label the rows of the heatmap by numbers
key<-data.frame(module=rownames(moduleTraitCor),id=seq(1:length(rownames(moduleTraitCor))),stringsAsFactors=FALSE)
# write.csv(key,file="ModuleNameColor_key_ModuleNumber.csv")

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
   ####
  
  moduleTraitCor2<-moduleTraitCor
  key.id2<-match(rownames(moduleTraitCor2),key$module)
  rownames(moduleTraitCor2)<-key$id[key.id2]
  
  
 # Display the correlation values within a heatmap plot
  
  if(nrow(moduleTraitCor2)>6){
   x.pv<-pvclust(moduleTraitCor2,nboot=100)
   colnames(moduleTraitCor2)<-gsub("Repetitive element","Rptv. Element",colnames(moduleTraitCor2))
  colnames(moduleTraitCor2)<-gsub("Endogenous Retrovirus","Endg. Retrovirus",colnames(moduleTraitCor2))
   #omit grey module.
  moduleTraitCor2<-moduleTraitCor2[which(key$module!="MEgrey"),]
   heatCor<-Heatmap(moduleTraitCor2,
                   column_names_gp=gpar(fontsize=10),
                   cluster_columns=x.pv$hclust,
                   cluster_rows=FALSE,
                   row_names_side="left",
                   name="correlation(x)",
                heatmap_legend_param=list(color_bar="continuous"),
                 column_title_gp=gpar(fontsize=13),
                 row_names_gp=gpar(fontsize=12),
                 column_title = expression("Associative Relationships of Gene Modules With Transposable Element Family Types  (*p.val" <= "0.05)"),
        cell_fun=function(j,i,x,y,w,h,col){
        weighted.Pvalue<-corPvalueStudent(moduleTraitCor2,nSamples)
                   if(weighted.Pvalue[i,j]<0.05){
                 #  grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                    grid.text("*",x,y-0.005 )
                    }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
           }
         )
     }else{
 heatCor<-Heatmap(moduleTraitCor,cluster_rows=FALSE,row_names_side="left",name="cor(x)", column_title = paste0("Module-Repeat ",how," Biotype relationships"))
  print(heatCor)
  }
  

#####add Heatcor of row annotation subset. 
  grey.id<-which(key$module=="MEgrey")
  toGo<-data.frame(module=unique(unlist(xN)),stringsAsFactors=FALSE)
  mC<-moduleTraitCor[rownames(moduleTraitCor)%in% paste0("ME",toGo$module),]
  mc.id<-match(rownames(mC),key$module)
  rownames(mC)<-key$id[mc.id]
 colnames(mC)<-gsub("Repetitive element","Rptv. Element",colnames(mC))
  colnames(mC)<-gsub("Endogenous Retrovirus","Endg. Retrovirus",colnames(mC))
  write.csv(mC,file="Module-Repeat-CPM-Biotype.relationships.csv")
  weighted.Pvalue<-corPvalueStudent(mC,nSamples)
 mC<- mC[which(rownames(mC)!=grey.id),]

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
  

###add activation direction
  activation.direction<-matrix(data=0,nrow=(nrow(moduleTraitCor2)+1),ncol=length(names(xN)))
  rownames(activation.direction)<-rownames(key)
  colnames(activation.direction)<-names(xN)
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

 activation.df<-sapply(df,function(x) as.numeric(x[grep("Total",rownames(x)),]$logFC))
  names(activation.df)<-paste0("ME",names(activation.df))
 dir.id<-match(names(activation.df),key$module)
  names(activation.df)<-dir.id
 dir.id2<-match(names(activation.df),rownames(activation.direction))

 activation.direction[dir.id2, which(colnames(activation.direction)==pathwaysToPick[i])]<-activation.df
}

########i
   activation.direction<-activation.direction[which(key$module!="MEgrey"),]

  activation.Heat<-Heatmap((activation.direction),col=colorRamp2(c(-2,0,2),c("black","white","orange")),name="Activity Level",column_title="Inferred Function",column_title_gp=gpar(fontsize=11),heatmap_legend_param=list(color_bar="continuous"),cluster_columns=FALSE,cluster_row=FALSE,show_row_names=FALSE, column_title_side="top" )

   par(mar = c(6, 14, 3, 3));
 # Display the correlation values within a heatmap plot
  draw(activation.Heat+heatCor)
  readkey()
#######
   ###fisher test mainly designed for 2 enrichment functions Immune and Inflammation not general
 pdf(paste0("Heat_correlation_Fisher",how,"_",p.value,"_plots.pdf"),width=12,height=12)
    par(mar = c(9, 10, 3, 3));
     print(activation.Heat+heatCor)
   dev.off()



  chi.df<-matrix(ncol=ncol(moduleTraitCor2),0)
  colnames(chi.df)<-colnames(moduleTraitCor2)
  chi.df<-as.data.frame(chi.df)
  t<-matrix(nrow=0,ncol=4,0)
  t2<-matrix(nrow=2,ncol=0,0)


  for(i in 1:ncol(moduleTraitCor2)){
     immune.up.Cor_pos<-sum(table(which(mC[names(which(activation.direction[,1]>0)),i]>0)))
     immune.up.Cor_neg<-sum(table(which(mC[names(which(activation.direction[,1]>0)),i]<0)))

     immune.down.Cor_pos<-sum(table(which(mC[names(which(activation.direction[,1]<0)),i]>0)))
     immune.down.Cor_neg<-sum(table(which(mC[names(which(activation.direction[,1]<0)),i]<0)))

     inf.up.Cor_pos<-sum(table(which(mC[names(which(activation.direction[,2]>0)),i]>0)))
     inf.up.Cor_neg<-sum(table(which(mC[names(which(activation.direction[,2]>0)),i]<0)))

     inf.down.Cor_pos<-sum(table(which(mC[names(which(activation.direction[,2]<0)),i]>0)))
     inf.down.Cor_neg<-sum(table(which(mC[names(which(activation.direction[,2]<0)),i]<0)))
    d<-rbind(c(immune.up.Cor_pos,immune.down.Cor_pos,inf.up.Cor_pos,inf.down.Cor_pos),c(immune.up.Cor_neg,immune.down.Cor_neg,inf.up.Cor_neg,inf.down.Cor_neg))
colnames(d)<-c("Immune.up","Immune.down","Inflam.up","Inflam.down")

    colnames(t)<-colnames(d)
    d<-data.frame(d)
    t<-rbind(t,d)
    rownames(t)[which(rownames(t)==1)]<-c(paste0(colnames(moduleTraitCor2)[i],".CorPos"))
    rownames(t)[which(rownames(t)==2)]<-c(paste0(colnames(moduleTraitCor2)[i],".CorNeg"))
    chi.df[,i]<-fisher.test(d)$p.value

    t2<-cbind(t2,d)
    colnames(t2)[which(colnames(t2)=="immune.pos")]<-paste0(colnames(moduleTraitCor2)[i],".",colnames(d)[1])
    colnames(t2)[which(colnames(t2)=="immune.neg")]<-paste0(colnames(moduleTraitCor2)[i],".",colnames(d)[2])
    colnames(t2)[which(colnames(t2)=="inf.pos")]<-paste0(colnames(moduleTraitCor2)[i],".",colnames(d)[3])
    colnames(t2)[which(colnames(t2)=="inf.neg")]<-paste0(colnames(moduleTraitCor2)[i],".",colnames(d)[4])
  }
    rownames(t2)<-c("positive.enriched","negative.enriched")
  global.correlation.chi<-chisq.test(colSums(t))$p.value
  pos<-t[grep(".CorPos",rownames(t)),]
   neg<-t[grep(".CorNeg",rownames(t)),]
  df.pos.neg<-data.frame(Pos=colSums(pos),Neg=colSums(neg))
  df.pos.neg<-t(df.pos.neg)
  pos.neg.correlation.chi<-chisq.test(df.pos.neg)$p.value
write.csv(t,file="Significant-Association.table.pairwise.Fisher.test.correlation.enrichment.csv")
  write.csv(chi.df,file="Fisher.pairwise.test.correlation.enrichment.csv")
  write.csv(global.correlation.chi,file="Global.chi.square.test.correlation.enrichment.csv")
  write.csv(pos.neg.correlation.chi,file="Positive.Negative.correlation.factor.correlation.enrichment.csv")
 # draw(activation.Heat+heatCor)
 # readkey()
 outputs<-t(chi.df)
 outputs<- cbind(outputs,p.adjust(outputs[,1],method="fdr",n=nrow(outputs)))
 
 outputs<-rbind(outputs,global.correlation.chi)
 outputs<-rbind(outputs,pos.neg.correlation.chi)
 colnames(outputs)<-c("p.value","FDR-adjust.p.value")
  return(outputs)


} ###main



