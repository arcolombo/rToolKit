#' @title renames repeat module colors to the top three correlations with biotypes
#' @description the rnames object is assigned color names, and we want to rename the repeat module colors after the top 3 highest correlated tx_biotype groups to assign an identity.  this will systematically rename the repeat module colors to leading biotype
#' @param rnames the object returned from wrcna method
#' @param return returns rnames object with traitCorRenamed object
#' @export
renameRepeatModuleColors<-function(rnames=NULL,rdbName="wrcnaDbLite.cpm.6.sqlite"){

  rbwModuleColors<-rnames[["moduleColors"]]
  rMEs<-rnames[["MEs"]]
  repeatModules<- listModuleColors(wgcnaDbLite(rdbName))
  rdatExpr<-rnames[["datExpr"]]
  rdatTraits<-rnames[["datTraits"]]
  annot<-rnames[["annot"]]
  stopifnot(rnames[["byWhich"]]=="repeat")
  allcolors<-dbListTables(dbconn(wgcnaDbLite(wgcnaDbName)))
  geneModules<-match.arg(geneModules,allcolors)
  repeatColors<-dbListTables(dbconn(wgcnaDbLite(wrcnaDbName)))
  stopifnot(all(repeatModules%in%repeatColors)==TRUE)
  key.id<-which(paste0("ME",geneModules)==colnames(MEs))
  moduleTraitCor<-rnames[["moduleTraitCor"]]
  
  rep.id<-match(rownames(rnames[["moduleTraitCor"]]),paste0("ME",repeatModules))
  ##grabs the top three leading correlation r.modules ~ tx-biotypes
  leading.terms<- lapply(rep.id,function(x) sort(rnames[["moduleTraitCor"]][x,],decreasing=TRUE)[1:3]  )
   names(leading.terms)<-rownames(rnames[["moduleTraitCor"]][rep.id,])
   newNames<-lapply(leading.terms,function(x) paste(names(x),collapse="-")  )
  repeat.key<-data.frame(leading.terms=unlist(newNames))
  repeat.key$leading.terms<-gsub(" ","_",repeat.key$leading.terms)
  write.csv(repeat.key,file="repeatModule.key.renamedLeadingBiotypes.correlations.csv")
  ##now to rename the correct object with the key
  #FIX ME: not sure which object to rename.  wish to rename the object used in moduleWiseAnalysis so that the output of module wise analysis prints out
  ##gene module - row number that matches the wgcna_Heatcor
  ##repeat module the repeat module will NOT print out the repeat color, but the leading 3 terms also in the wgcna_Heatcor summary (family) plot. 
 print("done.")
 return(repeat.key) 
}
