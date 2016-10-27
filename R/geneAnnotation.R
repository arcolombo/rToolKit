#' @title creates a geneInfo0 data frame which has gene symbol and module color
#' @description used mainly for WGCNA which returns a data frame of geneSym locus link ID module colr and correlations
#' @import biomaRt
#' @import arkas
#' @export
#' @return geneInof0 data frame
geneAnnotation<-function(datExpr,datTraits,species=c("Homo.sapiens","Mus.musculus")){

  stopifnot(any(grepl("ERCC",names(datExpr))==TRUE)==FALSE)

   if(class(datExpr)=="matrix"){
  datExpr<-as.data.frame(datExpr,stringsAsFactors=FALSE)
  }
  resValues<-names(datExpr)


require(biomaRt)
 species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
 if(species=="Homo.sapiens"){
  
  commonNomen<-"human"
  speciesMart<-.findMart(commonNomen)
  speciesSymbol<-"hgnc_symbol"  #hugo nomenclature human only 
         message("finding entrez IDs of top ensembl genes...")
         convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)

} else {
  
   commonNomen<-"mouse"
   speciesMart<-.findMart(commonNomen)
   speciesSymbol<-"mgi_symbol"
        message("finding entrez IDs of top ensembl genes...")
        convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
 }


  return(convertedEntrezID)
}


.findMart <- function(commonName=c("human","mouse"),host="www.ensembl.org"){#{{{

  dataset <- switch(match.arg(commonName),
                    human="hsapiens_gene_ensembl",
                    mouse="mmusculus_gene_ensembl")
  useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset, host=host)

} #}}}

