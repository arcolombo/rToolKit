#' @title runs qusage for pairwise or global runs
#' @description this is a point-wise qusage call not connected to database formations but used for kexp alone. Requires that hte comparison and controls be the leading column name able to be split by kexp2Group.  the qusageTables script is used to make a database. this script is used to make specific plots of qusage calls that target specific pathways.
#' @param p.cutoff threshold for FPR filtering alpha significance level
#' @import qusage
#' @import arkas
#' @import edgeR
#' @export
qusageResultsWithCI<-function(kexp,geneSetPath="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),how=c("cpm","tpm"),species=c("Homo.sapiens","Mus.musculus"),comparison=NULL,controls=NULL, showPlots=FALSE,comparisonNumber=1,pathwayTitle=NULL,paired=FALSE,read.cutoff=2,p.cutoff=0.05,asBand=TRUE,includeNames=FALSE,includeDE=FALSE){
 ###the plots picked here must match the pathways picked exactly by heat Cor
 #########
  ###TO DO: print everything to a table

  ##only controlling for specific comparison types for now  ###FIX ME: generalize
 geneSet<-MsigDB
  Pathwaytitle<-pathwayTitle
  #comparison<-match.arg(comparison,c("pHSC","Blast","RAEB","MDS"))
  comparison<-comparison 
  ##for comparisons with 2 groups, trios must exist
  if(comparisonNumber==2){
  kexp1<-kexpByTrio(kexp)
  phsc.lsc<-kexp2Group(kexp1,comparison=comparison,control=controls)
  }else{ ##if only 1 comparison group
   phsc.lsc<-kexp2Group(kexp,comparison=comparison,control=controls)
   }

  counts<-collapseBundles(phsc.lsc,"gene_name",read.cutoff=read.cutoff)
  ##BSN
  counts<-counts[,which(colSums(counts)>1)]
   dge<-DGEList(counts=counts)
   dge<-calcNormFactors(dge)
   expr<-cpm(dge,normalized.lib.sizes=TRUE,log=FALSE)
  expr2<-log2(1+expr)
  rownames(expr2)<-toupper(rownames(expr2))
  cN<-colnames(expr2)
  cN<-strsplit(cN,"_")
  labels<-unlist(lapply(cN,function(x) x[1]))

   contrast<-paste0(comparison,"-",controls)
   geneSets<-read.gmt(paste0(geneSetPath,geneSet))
  ##alphabetical order the geneSets must be
    geneSets<-geneSets[sort(names(geneSets))]

   if(paired==TRUE){
   pairs<-unlist(lapply(cN,function(x) x[2]))
   pairs.id<-match(toupper(pairs),toupper(pairs))
   qs.pHSC.results<-qusage(expr2,labels,contrast,geneSets,pairVector=pairs.id)
  }else if(paired==FALSE){
   qs.pHSC.results<-qusage(expr2,labels,contrast,geneSets)
  }
   max.mean<-max(qs.pHSC.results$path.mean)

  ##pathway indices : inflam  903, immune 394,410 ,wound 831, mapK 865,494, 
 
  qstab<-qsTable(qs.pHSC.results)
  t<-data.frame(qstab[,2:4])
  rownames(t)<-as.character(qstab$pathway.name)
    t<- t[order(rownames(t)),]
   alpha<-which(t$p.Value<=p.cutoff)
  rownames(t)[ which(t$p.Value<=p.cutoff)]<-paste0(rownames(t)[ which(t$p.Value<=p.cutoff)],"**")
    pchLabels<-as.numeric(rep(0,nrow(t)))
  rownames(t)<-tolower(rownames(t))
  rownames(t)<-paste0(toupper(substring(rownames(t),1,1)),substring(rownames(t),2))
  write.csv(t,file=paste0("qusage_",contrast,"_",MsigDB,"_enrichmentActivation.csv"))
  max.mean<-max(qs.pHSC.results$path.mean)
  par(cex.main=0.95,family='Helvetica') 
  plotDensityCurves(qs.pHSC.results,col=1:nrow(t), main=Pathwaytitle,xlab="Differential Enrichment Fold Change Level (Log2)",ylab="Distribution of Canonical Gene Set Activity" )
 legend("topleft",legend=rownames(t),col=1:nrow(t),pch=pchLabels,cex=0.99,bty='n',title=paste0("Significant Enrichment *p.val","\u2264",p.cutoff))
 mtext(paste0("Pairwise Comparison ",contrast),cex=0.95)
 # title(paste0("Differential Enrichment Activity of Apoptotic, Inflammatory, and Immune Canonical Gene Sets Comparing",comparison,"-",controls))
  
##Plot CIs
 x<-qsTable(qs.pHSC.results,number=1000)

 if(showPlots==TRUE){ 
  for(j in 1:nrow(t)){
 pathIndex<-j
  par(mar=c(10,4,2,2),mfrow=c(1,1),cex.main=0.98,family='Helvetica',cex.lab=0.9)
  x11(width=12,height=6)
  if(includeNames==FALSE){
   xLabels=NA
   }else{
   xLabels=NULL
  }
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast),pch=18,cex.xaxis=3.4,cex.axis=2,ylab="Differential Canonical Gene Set Activity (Log2 FC)",addGrid=FALSE,asBand=asBand,x.labels=xLabels )
#  readkey() 
  hea<-data.frame(mean=qs.pHSC.results$mean[qs.pHSC.results$pathways[[pathIndex]]][order(qs.pHSC.results$mean[qs.pHSC.results$pathways[[pathIndex]]],decreasing=TRUE)],SD=qs.pHSC.results$SD[qs.pHSC.results$pathways[[pathIndex]]][order(qs.pHSC.results$SD[qs.pHSC.results$pathways[[pathIndex]]],decreasing=TRUE)])
 write.csv(hea ,file=paste0("qusage.WelchTest.GeneSet.",contrast,".",names(qs.pHSC.results$pathways)[j],".csv"))
 ##its own page
###
  readkey()
   if(species=="Mus.musculus"){
   ##arkas annotation for mouse has gene symbols in lower case
   geneSets[[j]]<-firstup(tolower(geneSets[[j]]))
   }  
   imm<-phsc.lsc[rowRanges(phsc.lsc)$gene_name%in%geneSets[[j]],]
  gwa.imm<-geneWiseAnalysis(imm,design=metadata(imm)$design,how="cpm",adjustBy="BH",species=species)
   tpm<-collapseTpm(imm,"gene_name")
   path.Heat<-Heatmap(log(1+tpm[rownames(tpm)%in%as.character(gwa.imm$limmaWithMeta$Gene.symbol),] ),name="log(1+tpm)",column_title=paste0(names(geneSets)[j]," DE Genes Comparing ",contrast),row_names_gp=gpar(fontsize=9))
    write.csv(gwa.imm$limmaWithMeta,file=paste0(names(geneSets)[j],"_DEGenes_",contrast,".csv"))
  pdf(paste0(names(geneSets)[j],"_DEGenes_Heatmap_",contrast,".pdf"),height=14)
  print(path.Heat)
  dev.off()
  } ##individual path
  }##show plot

############PDF PLOT#################################################
 cairo_pdf(paste0("Qusage_Plotting_Results_Activation_",contrast,"_",MsigDB,".pdf"),family='Helvetica')
  par(cex.main=0.95)
  plotDensityCurves(qs.pHSC.results,col=1:nrow(t), main=Pathwaytitle,xlab="Differential Pathway Enrichment Level (Log2)",ylab="Distribution of Canonical Gene Set Activity",xlim=c(-max.mean-2,max.mean+2)  )
 legend("topleft",legend=rownames(t),col=c(1:nrow(t)),pch=pchLabels,bty='n',title=paste0("Significant Enrichment *p.val","\u2264",p.cutoff)  )
 mtext(paste0("Pairwise Comparison ",contrast),cex=0.95)
 dev.off()

 pdf(paste0("Qusage_Plotting_Results_CI_",contrast,"_",MsigDB,".pdf"),width=12,height=7)

###
 for(j in 1:nrow(t)){
 pathIndex<-j
  par(mar=c(12,5.3,2,2),mfrow=c(1,1),cex.main=1.2,family='Helvetica')
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast),pch=18,cex.xaxis=2.6,asBand=TRUE,col=2,ylab="Differential Canonical Gene Set Activity (Log2 FC)",x.labels=NA,addGrid=FALSE )

   ##its own page

   } ##individual path
 
  dev.off()
####################PDF PLOTTING##################################
 print("done.\n")
 return(x)
}#main

firstup <- function(x) {
   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
x
}
