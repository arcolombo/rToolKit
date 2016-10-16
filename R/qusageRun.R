#' @title runs qusage
#' @description a general method for calling qusage
#' @import qusage
#' @export
#' @return enrichment data set
qusageRun<-function(cnts_mt=NULL,gmt.path="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),comparison="pHSC",control="LSC",module=NULL,plotTable=FALSE,paired=TRUE){
 
##create labels from colnames
  cN<-colnames(cnts_mt)
  cN<-strsplit(cN,"_")
  labels<-unlist(lapply(cN,function(x) x[1]))
  contrast<-as.character(paste0(comparison,"-",control))


##create patient pairs 
   if(paired==TRUE){
   pairs<-unlist(lapply(cN,function(x) x[2]))
   pairs.id<-match(toupper(pairs),toupper(pairs))
  } 
##read MSigDB selection from input
  geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))
  geneSets<-read.gmt(paste0(gmt.path,"/",geneSet))
  cnts_pL2<-cnts_mt[which(rownames(cnts_mt)!=""),]
  ##collapse duplicated gene names with unique EnsIDs
    if(any(duplicated(rownames(cnts_pL2)))){
     stop("found duplicated hgnc names with differen EnsG ids please getGenes...\n")
      } 
 cnts_pL3<-cnts_pL2
  if(paired==TRUE){
 qs.results<-qusage(cnts_pL3,labels,contrast,geneSets,pairVector=pairs.id)
  } else{
  qs.results<-qusage(cnts_pL3,labels,contrast,geneSets)
  }
  #print(contrast)
  #print(qsTable(qs.results,number=40 ))
  if(plotTable==TRUE){
  par(mar = c(10, 6, 3, 3));
  plot(qs.results)
   title(paste0(contrast,".",module))
  pdf(paste0("qs_",module,"_",geneSet,"_",contrast,".pdf"))
  plot(qs.results)  
  dev.off()
  }
  return(qsTable(qs.results,number=40))
} ##main
