#' @title this will tabulate the enriched pathways as a proportion
#' @description the null hypothesis is that within each repeat biotype the proportions of inflammation response are identical to inflammation respsonse and non-inflammation response pathways the proportions are identical.  do all repeat biotypes associate equally to inflammation response?  this scripte will query the enrichment database , and tabulate the sum of inflammation related and non-inflammation releated enrichments.  note taht these enrichments were calculated as the most significantly cross-correlated genes for each partition module.  this is for an independent Chi-Square Test
#' @param qusageDb the dbname of the qusageDbLite
#' @param inflammatory the pathways related to MSigDb call
#' @param tx.Biotypes the cross correlated traits of interest
#' @param contrast phsc, or blast,  when tabulating phsc and blast data the test because paired,  when looking at individual contrasts we have independence
#' @export
#' @return a data frame of proportions between (inflammatory response,non inflammatory response) and proportions for each tx.Biotypes
repeatOdds<-function(qusageDb="qusageDbLite.cpm.sqlite",inflammatory=c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",  "HALLMARK_INTERFERON_GAMMA_RESPONSE",  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL2_STAT5_SIGNALING"),tx.Biotypes=c("Alu","DNA transposon","Endogenous Retrovirus","ERV1","ERV3","ERVK","ERVL","L1","L2","LTR Retrotransposon","Satellite")  ,contrast=c("phsc","blast"), paired=FALSE ){

  contrast<-match.arg(contrast,c("phsc","blast"))

  allcolors<-dbListTables(dbconn(qusageDbLite(qusageDb)))
  allcolors<-allcolors[!grepl("metadata",allcolors)]
  allcolors<-allcolors[!grepl("kexp",allcolors)]
    
  prop.df<-data.frame(inflammation_response=rep(0,length(tx.Biotypes)),
                    non_inflammation_response=rep(0,length(tx.Biotypes)))
  prop.df<-t(prop.df)
  colnames(prop.df)<-tx.Biotypes

  if(paired==FALSE){
##start tx bio for loop
  for(j in 1:length(tx.Biotypes)){
  prop.id<-grep(tx.Biotypes[j],colnames(prop.df))
  ##start color loop
    for(i in 1:length(allcolors)){
  enrichPath<-pathways(qusageDbLite(qusageDb),Module.color=allcolors[i],tx.Biotype=tx.Biotypes[j],contrast=contrast)
  enrichPath<-enrichPath[!is.na(enrichPath$FDR),]
  total<-nrow(enrichPath)
  infl.proportion<-nrow(enrichPath[enrichPath$pathway_name%in%inflammatory,])
  non.infl.proportion<-nrow(enrichPath[!enrichPath$pathway_name%in%inflammatory,])
  stopifnot(infl.proportion+non.infl.proportion==total)
  prop.df[1,prop.id]<-prop.df[1,prop.id]+ infl.proportion
  prop.df[2,prop.id]<-prop.df[2,prop.id]+ non.infl.proportion
  } #color loop
 }#tx bio loop


  prop.df<-rbind(prop.df,colSums(prop.df))
  prop.df<-cbind(prop.df,rowSums(prop.df)) 
  rownames(prop.df)[3]<-paste0(contrast,".total")
  colnames(prop.df)[ncol(prop.df)]<-paste0(contrast,".total")
  return(prop.df)
  }
if(paired==TRUE){
  ##now run a paired match 
  for(i in 1:length(tx.Biotypes)){
  prop.id<-grep(tx.Biotypes[j],colnames(prop.df))
  ##start color loop
    for(j in 1:length(inflammatory)){ 
      for(k in 1:length(allcolors)){
   phsc.enrichPath<-pathways(qusageDbLite(qusageDb),Module.color=allcolors[k],tx.Biotype=tx.Biotypes[j],contrast="phsc")
   blast.enrichPath<-pathways(qusageDbLite(qusageDb),Module.color=allcolors[i],tx.Biotype=tx.Biotypes[j],contrast="blast")

  enrichPath<-enrichPath[!is.na(enrichPath$FDR),]
  total<-nrow(enrichPath)
  infl.proportion<-nrow(enrichPath[enrichPath$pathway_name%in%inflammatory,])
  non.infl.proportion<-nrow(enrichPath[!enrichPath$pathway_name%in%inflammatory,])
  stopifnot(infl.proportion+non.infl.proportion==total)
  prop.df[1,prop.id]<-prop.df[1,prop.id]+ infl.proportion
  prop.df[2,prop.id]<-prop.df[2,prop.id]+ non.infl.proportion
  } #color loop
 } #each inflam path
 }#tx bio loop


   } #paired  
} ##main 


matchPair<-function(){


}
