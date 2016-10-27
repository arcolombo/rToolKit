#' @title queries a data frame for ISGs,NK ligands and DDX helicases
#' @description this is a first draft to query a data frame that has a column df gene names to return information about ISG, NK ligands or DDX helicases 
#' @param daFr should either be a driver data frame of arkas GWA df
#' @param geneSymbol_col is the column named Gene.symbol from arkas of hgnc_symbol from repeatToolKit drivers 
#' @export
#' @return information only, no object
queryDF<-function(daFr=NULL,geneSymbol_col=NULL,NKLigs=c("CLEC2B","KIR3DL3","CLEC4F","CLEC2A","KLRG2","CLEC4G","ULBP3","RAET1E","RAET1L","CLEC4C","CLEC2L","KLRD1","CLEC10A","CLEC3A","CLEC3B","CLEC19A","CLEC1A","MICA","RAET1G"),ISG=c("IFNGR1","TLR2","IFNAR2","IFITM3","IFNL4","IFNLR1","TLR9","IFNL1","CCL16","TRIL","CCL4","TLR5","IFNL2","IFNL3","CCL17","CCL22","TLR7","TLR10","IFITM5","IFNB1","CCL25","CCL15","IFRD2","ISG20"),ISG.pre=c("IL","DDX","IFN","TLR","IFIT","CCL","ISG","IRF","NFKB"),NK.pre=c("CLEC","KIR","KLR","RAET","MIC","CLEX","ULBP")){


  if(is.null(geneSymbol_col)==TRUE){
   id<-grep("hgnc_symbol",colnames(daFr))
   if(length(id)==0){
    id<-grep("Gene.symbol",colnames(daFr))
    stopifnot(length(id)>0)
    }  
  }
  targets<-c(ISG.pre,NK.pre)
  ##Fix me print a report of query. 
  for(i in 1:length(targets)){
  seed<-daFr[grep(targets[i], as.character(daFr[,id])),]
   if(nrow(seed)>0){
   print(seed)  
  readkey()
  }
  }


 } #main 
