#' @title given the module of interest printed to a csv, we can read that data, and analyze the pathway behavior of specific pathways
#' @description reads the qusage module pathway activity and analyzes specific behvaviors within that function
#' @param db the wgcnaDbName
#' @param qdb the qusageDbLite name
#' @param functionKeyWords this is a single character that will be grep'd from the list of pathway modules 
#' @export
#' @return a data frame with the queried keyword in every module and pathway information
moduleFunctionDataFrame<-function(db=NULL,qdb=NULL,functionKeyWords=NULL){
 #task search through all the geneModules for key Words and aggregate the logFc into analysis.
  allcolors<-listModuleColors(wgcnaDbLite(db))
  df<-data.frame()
  for(colR in allcolors){
  t<-pickPathway(qusageDbLite(qdb),Module.color=colR,keyWord=functionKeyWords)
  if(any(is.na(t$pvalue))==FALSE){
  df<-rbind(df,t)
   }
  }
  df$weight<-df$logFC*df$keyWord_Percentage 
 return(df)

}
