#' @title gets gene name counts
#' @description this will input a data frame of ENSG counts and find the entrez and hgnc gene names and then collapse by gene name.  this is used for qusageTablesFromWGCNA to collapse by tmm normalized ensGID counts.  arkas collapesBunders(kexp,gene_name) will give you counts(..) of gene_name bundles which is not exactly the same as using tmm normalized ENSGIDs and then collapsing tmm ENSGID counts by gene name.  slight difference but qusage requires it.
#' @param df a tmm normalized ensGID data frame
#' @param species Homo sapiens for now
#' @export 
getGenes<-function(df=NULL,species="Homo.sapiens"){
  df<-as.data.frame(df)
  gxs<-data.frame(gene_id=rownames(df),stringsAsFactors=TRUE)
  rownames(gxs)<-as.character(gxs$gene_id)
  entrezID <- getEntrezIDs(gxs, species)
  entrez.match<-match(rownames(gxs),names(entrezID))
  stopifnot(all(rownames(gxs)==names(entrezID)[entrez.match]))
  gxs$entrezid<-entrezID[entrez.match]
  symbolNames<-getSymbols(gxs,species)
  symbol.match<-match(rownames(gxs),names(symbolNames))
  stopifnot(all(rownames(gxs)==names(symbolNames)[symbol.match]))
  gxs$hgnc_symbol<-symbolNames[symbol.match]
  id<-match(rownames(df),gxs$gene_id)
  df2<-data.frame(df,hgnc_symbol=as.factor(gxs[id,3]))
 stopifnot(rownames(df)==rownames(gxs[id,]))
 df3<-df2[!is.na(df2$hgnc_symbol),]
 df4<-split.data.frame(df3,df3$hgnc_symbol)
  totalCol<-ncol(df)
  bundled<-do.call(rbind,
               lapply(df4,function(x) colSums(x[,c(1:totalCol)]))
               )
  stopifnot(any(duplicated(rownames(bundled)))==FALSE)
  return(bundled)
}
