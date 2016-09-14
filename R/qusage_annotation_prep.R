#' @title creates a kexp data frame which has gene symbol rownames
#' @description used mainly for qusage calls which returns a data frame of gene counts with rownames as gene HUGO names
#' @import biomaRt
#' @import arkas
#' @export
#' @return gene name data frame
qusage_annotation_prep<-function(kexp,species=c("Homo.sapiens","Mus.musculus")){

  datExpr<-collapseBundles(kexp,"gene_id")
  resValues<-rownames(datExpr)


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

 ###annotate rowRanges with gene name then collapse by gene_name
   ######add gene names
  id<-which(convertedEntrezID$hgnc_symbol!="")
  conv<-convertedEntrezID[id,]
  id2<-match(conv$ensembl_gene_id,rowRanges(kexp)$gene_id)
  rowRanges(kexp)$gene_name[id2]<-convertedEntrezID$hgnc_symbol
 return(kexp)
} # main 


.findMart <- function(commonName=c("human","mouse"),host="www.ensembl.org"){#{{{

  dataset <- switch(match.arg(commonName),
                    human="hsapiens_gene_ensembl",
                    mouse="mmusculus_gene_ensembl")
  useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset, host=host)

} #}}}

