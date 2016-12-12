#' @title runs qusage for pairwise or global runs
#' @description this is a point-wise qusage call not connected to database formations but used for kexp alone. Requires that hte comparison and controls be the leading column name able to be split by kexp2Group.  the qusageTables script is used to make a database. this script is used to make specific plots of qusage calls that target specific pathways.
#' @import qusage
#' @import arkas
#' @import edgeR
#' @export
qusageResults<-function(kexp,geneSetPath="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),how=c("cpm","tpm"),species=c("Homo.sapiens"),comparison=c("pHSC","Blast"),controls="LSC", keyWords=c("immun","inflam","apopto","death","kappab","wound"),showPlots=FALSE,comparisonNumber=1,paired=FALSE,read.cutoff=2) {
 ###the plots picked here must match the pathways picked exactly by heat Cor
 ## pp<- pickPathway(qusageDbLite(qdbname),keyWord="immun")
  ##pp2<-sapply(pp,function(x) unique(x$pathway_name))
  ##unique_names <-unique(unlist(pp2))
  ##uni_names<-unique_names[!grepl(keyWords,unique_names)]
  ##uni_names has all the pathway activations for keyWords i that is used in calculating the activation direction in HeatCor image. 
  ##fix me: add a loop here to programmatically query the stage-level pathway analysis to match/search for the pathways found in the module pathway analysis
 #########
  ###TO DO: print everything to a table


 geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))
  ##only controlling for specific comparison types for now  ###FIX ME: generalize
 comparison<-match.arg(comparison,c("pHSC","Blast","RAEB","MDS"))
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

   if(paired==TRUE){
   pairs<-unlist(lapply(cN,function(x) x[2]))
   pairs.id<-match(toupper(pairs),toupper(pairs))
   qs.pHSC.results<-qusage(expr2,labels,contrast,geneSets,pairVector=pairs.id)
  }else if(paired==FALSE){
   qs.pHSC.results<-qusage(expr2,labels,contrast,geneSets)
  }
  ##pathway indices : inflam  903, immune 394,410 ,wound 831, mapK 865,494, 
  if(showPlots==TRUE){
   plot(qs.pHSC.results)
   title(paste0("Global Pathway ",comparison,"-",controls))
   readkey()
   }

 x<-qsTable(qs.pHSC.results,number=1000)
 if(showPlots==TRUE){
  for(i in keyWords  ){
  #in
  pathIndex<-as.numeric(rownames(x[grep(i,x$pathway.name,ignore.case=TRUE),]))
  if(length(pathIndex)==0){
   next
   } 
 sigPath<-x[grep(i,x$pathway.name,ignore.case=TRUE),]
   sigId<-which(sigPath$p.Value<=0.05) 
 tp3<-as.character((x[grep(i,x$pathway.name,ignore.case=TRUE),])$pathway.name)
  tp3<-gsub("OF_","",tp3)
  if(length(sigId)>0){
  tp3[sigId]<-paste0(tp3[sigId],"**")
  }
 if(i=="immun"){
  tp3<-gsub("IMMUNE_","",tp3)
  legendTitle<-"Immune Pathways ** pvalue<=0.05"
  }else if(i=="inflam"){
  legendTitle<-"Inflammatory Pathways ** pvalue<=0.05"
  }else if(i=="apopto"){
  legendTitle<-"Apoptotic Pathways ** pvalue<=0.05"
  }else if(i=="death"){
  legendTitle<-"Programmed Cell Death ** pvalue<=0.05"
  }else if(i=="kappab"){
  legendTitle<-"NF Kappab ** pvalue<=0.05"
  }else if(i=="wound"){
  legendTitle<-"Wound Healing ** pvalue<=0.05"
  }
 plotDensityCurves(qs.pHSC.results,path.index=pathIndex,xlim=c(-2,2),col=1:length(pathIndex),main=paste0("Enriched Immune System Pathway ",contrast)  )
 legend("topleft",legend=tp3,col=1:length(tp3),pch=0.8,cex=0.6,bty='n',title=legendTitle)
  readkey()
  for(j in pathIndex){
  par(mar=c(12,5.3,2,2),mfrow=c(1,2),cex.main=0.44)
  plotDensityCurves(qs.pHSC.results,path.index=j,xlim=c(-2,2),main=paste0(x[which(rownames(x)== j),1]," ",contrast))
  legend("topleft",legend=signif(x[which(rownames(x)==j),3],2),bty='n',cex=0.6,pch=0.8,title="p.value" )
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast) )
  readkey()
 ##its own page
 shortTitle<-unlist(strsplit(as.character(x[which(rownames(x)== j),1]),"_"))
 if(length(shortTitle)>=3){
  keyWord.Id<-!grepl("of",shortTitle,ignore.case=TRUE)
  shortTitle<-shortTitle[keyWord.Id]
  keyWord.Id2<-!grepl("by",shortTitle,ignore.case=TRUE)
  shortTitle<-shortTitle[keyWord.Id2] 
  if(length(shortTitle)>=4){
 shortTitle<-paste0(shortTitle[1],"_",shortTitle[2],"_",shortTitle[3],"_",shortTitle[4] )
  }else{
shortTitle<-paste0(shortTitle[1],"_",shortTitle[2],"_",shortTitle[3] )
  }
 }else{
  shortTitle<-paste0(shortTitle[1],"_",shortTitle[2])
 } 
 plotGeneSetDistributions(qs.pHSC.results,path.index=j ,groupLabel=contrast,main=shortTitle )
###
  readkey()
  } ##individual path
}## merged  
}##show plots

############PDF PLOT#################################################
 pdf(paste0("Qusage_Plotting_Results_",contrast,".pdf"))
   plot(qs.pHSC.results)
   title(paste0("Global Pathway ",contrast))
 for(i in keyWords  ){
  #in
  pathIndex<-as.numeric(rownames(x[grep(i,x$pathway.name,ignore.case=TRUE),]))
  if(length(pathIndex)==0){
  next
  }
  sigPath<-x[grep(i,x$pathway.name,ignore.case=TRUE),]
 write.csv( sigPath,file=paste0("QuSAGE_",i,"_",comparison,".vs.",controls,".csv"))
 
  sigId<-which(sigPath$p.Value<=0.05)
 tp3<-as.character((x[grep(i,x$pathway.name,ignore.case=TRUE),])$pathway.name)
  tp3<-gsub("OF_","",tp3)
  if(length(sigId)>0){
  tp3[sigId]<-paste0(tp3[sigId],"**")
  }
 if(i=="immun"){
  tp3<-gsub("IMMUNE_","",tp3)
  legendTitle<-"Immune Pathways ** pvalue<=0.05"
  }else if(i=="inflam"){
  legendTitle<-"Inflammatory Pathways ** pvalue<=0.05"
  }else if(i=="apopto"){
  legendTitle<-"Apoptotic Pathways ** pvalue<=0.05"
  }else if(i=="death"){
  legendTitle<-"Programmed Cell Death ** pvalue<=0.05"
  }else if(i=="kappab"){
  legendTitle<-"NF Kappab ** pvalue<=0.05"
  }else if(i=="wound"){
  legendTitle<-"Wound Healing ** pvalue<=0.05"
  }
 plotDensityCurves(qs.pHSC.results,path.index=pathIndex,xlim=c(-2,2),col=1:length(pathIndex),main=paste0("Enriched Immune System Pathway ",contrast)  )
 legend("topleft",legend=tp3,col=1:length(tp3),pch=0.8,cex=0.6,bty='n',title=legendTitle)
 
###
 for(j in pathIndex){
  par(mar=c(12,5.3,2,2),mfrow=c(1,2),cex.main=0.44)
  plotDensityCurves(qs.pHSC.results,path.index=j,xlim=c(-2,2),main=paste0(x[which(rownames(x)== j),1]," ",contrast))
  legend("topleft",legend=signif(x[which(rownames(x)==j),3],2),bty='n',cex=0.6,pch=0.8,title="p.value" )
  plotCIsGenes(qs.pHSC.results,path.index=j, main=paste0(x[which(rownames(x)== j),1]," ",contrast) )
   ##its own page
 shortTitle<-unlist(strsplit(as.character(x[which(rownames(x)== j),1]),"_"))
 if(length(shortTitle)>=3){
  keyWord.Id<-!grepl("of",shortTitle,ignore.case=TRUE)
  shortTitle<-shortTitle[keyWord.Id]
  keyWord.Id2<-!grepl("by",shortTitle,ignore.case=TRUE)
  shortTitle<-shortTitle[keyWord.Id2]
  if(length(shortTitle)>=4){
 shortTitle<-paste0(shortTitle[1],"_",shortTitle[2],"_",shortTitle[3],"_",shortTitle[4] )
  }else{
shortTitle<-paste0(shortTitle[1],"_",shortTitle[2],"_",shortTitle[3] )
  }
 }else{
  shortTitle<-paste0(shortTitle[1],"_",shortTitle[2])
 }
 plotGeneSetDistributions(qs.pHSC.results,path.index=j ,groupLabel=contrast,main=shortTitle )
###
   } ##individual path
 } ##merged 
  dev.off()
####################PDF PLOTTING##################################
 print("done.\n")
 save(qs.pHSC.results,file=paste0(comparison,".vs.",controls,".qusageResults_CPM.RData"),compress=TRUE)
 return(x)

}#main
