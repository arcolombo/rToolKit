#' @title returns the central hub genes for analysis within a kexp
#' @description performs a connectivity analysis and returns the top connected genes within an experiment.
#' @param kexp kexp that is annotated
#' @param lnames wgcna.R output
#' @param topNumber  connectivity cutoff
#' @import WGCNA
#' @import arkas
#' @export
#' @return a list of genes and their connectivity
topHubGenes<-function(kexp,lnames=NULL,topNumber=20){
  datExpr1g<-lnames[["datExpr"]]
  gm1<-signedKME(datExpr1g,lnames[["MEs"]])
  colnames(gm1)<-gsub("kME","",colnames(gm1))
  MMPvalue1<-corPvalueStudent(as.matrix(gm1),dim(datExpr1g)[[2]])
  colnames(MMPvalue1)<-paste0("p.val.",colnames(MMPvalue1))
modulesA1<-lnames[["moduleColors"]]
    Gene<-colnames(datExpr1g)
   kMEtable1<-cbind(Gene,modulesA1)
   colnames(kMEtable1)<-c("GeneID","Module")
  kMEtable1<-cbind(kMEtable1,gm1)
  all(kMEtable1[,1]==rownames(gm1))
  all(rownames(kMEtable1)==rownames(MMPvalue1))
  kMEtable1<-cbind(kMEtable1,MMPvalue1)
 ##now the gene modules comparison are to be compared with a reference using the assignments in the comparison
  

  #in order to compare gene module membership across experiments the ME must have the same samples across experiments.

#### within gene module top rank calculations
  topGenesKME<-NULL
  for(i in 1:ncol(gm1)){
  if(i%%31==0){
    cat(paste(".",""))
  }
  gm.id<-which(colnames(kMEtable1)==colnames(gm1)[i])
  gm.pvalue.id<-which(colnames(kMEtable1)==paste0("p.val.",colnames(gm1)[i]))
  kMErank1<-rank(gm1[,i],ties.method="average")
  index.range<-NULL
   ##pick the 1:topNumber of top connected gene members in a module. 
  ##unique ranks will filter out the ties.
  for(j in 1:length(kMErank1)){
  rank.range<-unique(sort(kMErank1,decreasing=TRUE))[j]
   ##we use unique to find the range of the ranks
   ##if NA
   if(j%%71==0){
    cat(paste(".",""))
  }

   if(is.na(rank.range)==TRUE){
    next
   }
  if(length(index.range)>topNumber){
   break
   }
   #now we find which kMErank1 elements match the jth rank in the range
      ind.x<-which(kMErank1==rank.range) 
     for(k in 1:length(ind.x)){
  ## ensure that the top ranked members were actually assigned to the module.
      if(kMEtable1[ind.x[k],]$Module==colnames(gm1)[i]){
          index.range<-cbind(index.range,ind.x[k])
          }
      } #for all indices found
    } #for j topNumber ranked
   ##find gene names
  tt<- kMEtable1[as.vector(index.range),c(gm.id,gm.pvalue.id)]
  t<-cbind(tt,rowRanges(kexp)[match(rownames(tt),rowRanges(kexp)$gene_id)]$gene_name)
  colnames(t)[3]<-"Gene_Name"
  topGenesKME[[i]]<-t
  names(topGenesKME)[i]<-colnames(gm1)[i]
  } #for all gm1
  lapply(topGenesKME,function(x) write.table(x,file="topHubGenes.AML.csv",append=TRUE,sep=',')) 

###########
 return(topGenesKME)
}
