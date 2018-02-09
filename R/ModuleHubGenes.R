#' @title returns the central hub genes for analysis within a kexp for a selected module color. topHubGenes finds hubs across all modules.
#' @description performs a connectivity analysis and returns the top connected genes within an experiment.
#' @param kexp kexp that is annotated
#' @param lnames wgcna.R output
#' @param topNumber  connectivity cutoff
#' @import WGCNA
#' @import arkas
#' @export
#' @return a list of genes and their connectivity
ModuleHubGenes<-function(kexp,lnames=NULL,module.color=NULL,connectivity.cutoff=0.6){
  stopifnot(is.null(module.color)==FALSE)
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
 
  df<-data.frame(kMEtable1[,c(which(colnames(kMEtable1)==module.color),which(colnames(kMEtable1)==paste0("p.val.",module.color)))])  
  t.hub<-df[which(df[,1]>connectivity.cutoff),]
blast<-kexp
   t.DF<-data.frame(id=rowRanges(blast[rowRanges(blast)$gene_id%in%rownames(t.hub),])$gene_id)
  t.DF2<- c(rowRanges(blast[rowRanges(blast)$gene_id%in%rownames(t.hub),])$gene_name)
  t.DF3<-(cbind(t.DF,t.DF2))
  id<-match(rownames(t.hub),as.character(t.DF3$id))
   t.DF4<-(t.DF3[id,])
   t.DF5<-(cbind(t.DF4,t.hub))
  t.DF6<-t.DF5[order(t.DF5[,3],decreasing=T),]
 ##now the gene modules comparison are to be compared with a reference using the assignments in the comparison
  

  #in order to compare gene module membership across experiments the ME must have the same samples across experiments.

###########
 return(t.DF6)
}
