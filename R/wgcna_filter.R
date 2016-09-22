#' @title filters low mean counts 
#' @description often it is needed to filter counts, this does that
#' @param datExpr from the wgcna object lnames
#' @param annot from the gene annotation call
#' @export
#' @import WGCNA
#' @return a data frame with p adjusted q value genes
wgcna_filter<-function(lnames,geneTraitModuleDF,datExpr,annot,qvalCut=0.1){


##the geneTraitModuleDF has pValues, there are not Fisher pvalues....  maybe add ?
##has a qvalue for each module color
q.list<-lapply(geneTraitModuleDF[,grepl("p.MM",colnames(geneTraitModuleDF))],
                function(x) qvalue(x)$qvalues)
names(q.list)<-sub("^p.","q.",names(q.list))

###module list p and q value
colors<-sub("^MM","",colnames(geneTraitModuleDF)[grepl("^MM",colnames(geneTraitModuleDF))])
###need to subset each individual module color genes,corrl,p.value,fisherp,qvalue  in a list.  filter and then write out to text files each module
moduleSig<-list()
  for(i in 1:length(colors)){
  modID<-which(colors[i]==substring(colnames(geneTraitModuleDF),3) )
  pID<-which(colors[i]==substring(colnames(geneTraitModuleDF),5))
  #print(colnames(geneTraitModuleDF)[c(modID,pID)])
  moduleSig[[i]]<-geneTraitModuleDF[,c(modID,pID)]
  qvalu<-qvalue(moduleSig[[i]][,2])$qvalues
  moduleSig[[i]]<-cbind(moduleSig[[i]],qvalu)
  names(moduleSig)[i]<-colors[i]
}
qFilter<-lapply(moduleSig,function(x)
               x[which(x$qvalu<qvalCut),])
#######################################
message("table under q.value ",qvalCut)
print(lapply(moduleSig,function(x) table(x[3]<qvalCut)))

 
###find the geneTrait pVal and qVal


 GeneName=colnames(datExpr)
  hugoID<-match(GeneName,annot$ensembl_gene_id)
  HugoName<-annot$hgnc_symbol[hugoID]
  StandardGeneScreeningResults=data.frame(GeneName,
                               PearsonCorrelation=GS1, 
                               p.Standard,
                               q.Standard,
                               HGNC=HugoName)

return(StandardGeneScreeningResults)
} 
