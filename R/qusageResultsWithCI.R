#' @title runs qusage for pairwise or global runs
#' @description this is a point-wise qusage call not connected to database formations but used for kexp alone. Requires that hte comparison and controls be the leading column name able to be split by kexp2Group.  the qusageTables script is used to make a database. this script is used to make specific plots of qusage calls that target specific pathways.
#' @param p.cutoff threshold for FPR filtering alpha significance level
#' @import qusage
#' @import arkas
#' @import edgeR
#' @export
qusageResultsWithCI<-function(kexp,geneSetPath="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),how=c("cpm","tpm"),species=c("Homo.sapiens"),comparison=NULL,controls=NULL, showPlots=FALSE,comparisonNumber=1,pathwayTitle=NULL,paired=FALSE,read.cutoff=2,p.cutoff=0.05){
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
  ##pathway indices : inflam  903, immune 394,410 ,wound 831, mapK 865,494, 
  if(showPlots==TRUE){
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
  plotDensityCurves(qs.pHSC.results,col=1:nrow(t), main=Pathwaytitle,xlab="Differential Enrichment Level",ylab="Distribution of Canonical Gene Set Activity",xlim=c(-max.mean-2,max.mean+2)  )
 legend("topleft",legend=rownames(t),col=1:nrow(t),pch=pchLabels,cex=0.99,bty='n',title=paste0("Significant Enrichment *p.val","\u2264",p.cutoff))
 mtext(paste0("Pairwise Comparison ",contrast),cex=0.95)
 # title(paste0("Differential Enrichment Activity of Apoptotic, Inflammatory, and Immune Canonical Gene Sets Comparing",comparison,"-",controls))
  
##Plot CIs
 x<-qsTable(qs.pHSC.results,number=1000)

  
  for(j in 1:nrow(t)){
 pathIndex<-j
  par(mar=c(12,5.3,2,2),mfrow=c(1,1),cex.main=0.98,family='Helvetica',cex.lab=0.9)
  x11(width=12,height=6)
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast),pch=18,cex.xaxis=1.2,ylab="Differential Canonical Gene Set Activity" )
#  readkey()
 write.csv(qs.pHSC.results$mean[qs.pHSC.results$pathways[[pathIndex]]][order(qs.pHSC.results$mean[qs.pHSC.results$pathways[[pathIndex]]],decreasing=TRUE)] ,file=paste0("qusage.GeneSet.",contrast,".",names(qs.pHSC.results$pathways)[j],".csv"))
 ##its own page
###
  readkey()
  } ##individual path
  }##show plots

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
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast),pch=18,cex.xaxis=2.6,asBand=TRUE,col=2,ylab="Differential Canonical Gene Set Activity" )

   ##its own page

   } ##individual path
 
  dev.off()
####################PDF PLOTTING##################################
 print("done.\n")
 return(x)
}#main
