#' @title filters low mean counts 
#' @description often it is needed to filter counts, this filters geneIDs from modules based on student-t p values.  we do have fisher pvalues, but not used to filter (yet)
#' @param lnames the object called from wgcna.R 
#' @param geneTraitModuleDF the object called from wgcna_scatterMod.R
#' @export
#' @import WGCNA
#' @return a data frame with p adjusted q value genes
wgcna_filter<-function(lnames,geneTraitModuleDF,qvalCut=0.1){
 message("found wgcna.dataInput")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  MEs<-lnames[["MEs"]]
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
  modulePvalFisher<-lnames[["modulePvalFisher"]]


##has a qvalue for each module color
q.list<-lapply(geneTraitModuleDF[,grepl("^p.MM",colnames(geneTraitModuleDF))],
                function(x) qvalue(x)$qvalues)
names(q.list)<-sub("^p.","q.",names(q.list))

###module list p and q value
colors<-sub("^MM","",colnames(geneTraitModuleDF)[grepl("^MM",colnames(geneTraitModuleDF))])
traits<-sub("^GCor.","",colnames(geneTraitModuleDF)[grepl("^GCor.",colnames(geneTraitModuleDF))])
### subsets each individual module color in terms of genes,corrl,p.value,fisherp,qvalue  in a list.  filter and then write out to text files each module
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
module.qFilter<-lapply(moduleSig,function(x)
               x[which(x$qvalu<qvalCut),])
#######################################
message("table under q.value ",qvalCut)
print(lapply(moduleSig,function(x) table(x$qvalu<qvalCut)))

##annotate qFilter
for(i in 1:length(module.qFilter)){
q.ID<-match(rownames(module.qFilter[[i]]),annot$ensembl_gene_id)
module.qFilter[[i]]<-cbind(module.qFilter[[i]],annot[q.ID,])
}
 
###find the geneTrait pVal and qVal
##this grabs student pvalues and fisher pvalues
traitSig<-list()
  for(i in 1:length(traits)){
  modID<-which(traits[i]==substring(colnames(geneTraitModuleDF),6) )
  pID<-which(traits[i]==substring(colnames(geneTraitModuleDF),8))
  #print(colnames(geneTraitModuleDF)[c(modID,pID)])
  traitSig[[i]]<-geneTraitModuleDF[,c(modID,pID)]
   ##
  studentPval.id<-grepl("p.GCor.",colnames(traitSig[[i]]))
  fisherPval.id<-grepl("pf.GCr.",colnames(traitSig[[i]]))
  student.qvalu<-qvalue(traitSig[[i]][,studentPval.id])$qvalues
  fisher.qvalu<-qvalue(traitSig[[i]][,fisherPval.id])$qvalues

  traitSig[[i]]<-cbind(traitSig[[i]],student.qvalu)
  traitSig[[i]]<-cbind(traitSig[[i]],fisher.qvalu)
  names(traitSig)[i]<-traits[i]
}
##filter based on student for now
trait.qFilter<-lapply(traitSig,function(x)
               x[which(x$student.qvalu<qvalCut),])

message("table trait significance under q.value ",qvalCut)
print(lapply(traitSig,function(x) table(x$student.qvalu<qvalCut)))

## annotate geneTrait 
##annotate qFilter
for(i in 1:length(trait.qFilter)){
q.ID<-match(rownames(trait.qFilter[[i]]),annot$ensembl_gene_id)
trait.qFilter[[i]]<-cbind(trait.qFilter[[i]],annot[q.ID,])
}

save(trait.qFilter,file="trait.qFilter.RData",compress=TRUE)
save(module.qFilter,file="module.qFilter.RData",compress=TRUE)
## print to csv?  some have 14,000 entries .... too many to print to csv


  StandardGeneScreeningResults=list(trait.qFilter=trait.qFilter,module.qFilter=module.qFilter)

return(StandardGeneScreeningResults)
} 
