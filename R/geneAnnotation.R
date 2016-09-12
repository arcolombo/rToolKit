#' @title creates a geneInfo0 data frame which has gene symbol and module color
#' @description used mainly for WGCNA which returns a data frame of geneSym locus link ID module colr and correlations
#' @import biomaRt
#' @export
#' @return geneInof0 data frame
geneAnnotation<-function(datExpr,datTraits,species="Mus.musculus"){
##FIX ME: locus link IDs could be added although EnsgID should be unambiguious

###FIX ME: entrez gene IDs must be added.
require(biomaRt)
  stopifnot(any(grepl("ERCC",names(datExpr))==TRUE)==FALSE)
 speciesMart<-useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
  if(class(datExpr)=="matrix"){
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  }
 convertedEntrezID<-getBM(filters="ensembl_gene_id",
 attributes=c("ensembl_gene_id","entrezgene",speciesSymbol="mgi_symbol",
 "description"),
 values=names(datExpr),
 mart=speciesMart)

return(convertedEntrezID)
}
