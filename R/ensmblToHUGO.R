#' @title map ensembl id to hugo via biomaRt
#' @description sometimes you just want gene names from a TPM list.  this uses a collapseTPM call to map to HUGO to allow for exploration
#' @import biomaRt
#' @param kexp a grimes kexp at the gene or repeat level
#' @param commonNomen character moues or human
#' @export
#' @return a list with the meta data added
ensmblToHUGO<-function(kexp,commonNomen=c("human","mouse"),byWhich=c("gene","repeat") ){
  commenNomen<-match.arg(commonNomen,c("human","mouse"))
   if(byWhich=="gene"){
   tpm<-collapseTpm(kexp,"gene_id")
   } else {
   kexp<-findRepeats(kexp)
   tpm<-collapseTpm(kexp,"tx_id")
    }
   
   resValues<-rownames(tpm)
 
 if (commonNomen=="human") {
   speciesMart<-.findMart(commonNomen)
    speciesSymbol<-"hgnc_symbol"  #hugo nomenclature human only 
         message("finding entrez IDs of top ensembl genes...")
         convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)
   }#human
  if(commonNomen=="mouse"){
   speciesMart<-.findMart(commonNomen)
   speciesSymbol<-"mgi_symbol"
        message("finding entrez IDs of top ensembl genes...")
        convertedEntrezID<-getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id","entrezgene",speciesSymbol,"description"),
                    values=resValues, #fitBundles ensembl Gene Ids
                    mart=speciesMart)


} #mouse
   converted<-convertedEntrezID
  
 mapped<-.formatWithMeta(tpm,converted,kexp)
 return(mapped)
} ##main
.findMart<-function(commonName=c("human","mouse"),host="www.ensembl.org"){#{{{

 dataset <- switch(match.arg(commonName),
                    human="hsapiens_gene_ensembl",
                    mouse="mmusculus_gene_ensembl")
  useMart("ENSEMBL_MART_ENSEMBL", dataset=dataset, host=host)
}






.formatWithMeta<-function(tpm,converted,kexp){ # {{{ format limma results
  res<-list()
  res$top<-tpm
  # create csv of limma counts, gene names, ensembl ID, biotypes; store into res
  index<-vector()
  for(i in 1:nrow(converted)){
    cols <- grep("ensembl_gene_id", colnames(converted))
    index[i] <- which(rownames(tpm) == converted[i, cols])
  } #indexing converted

  limmad <- res$top[index,]
  limmad <- cbind(limmad,
                  converted[, grep("entrezgene",colnames(converted))],
                  converted[, grep("_symbol",colnames(converted))],
                  converted[, grep("ensembl_gene_id",colnames(converted))],
                  converted[, grep("description",colnames(converted))] )
  colnames(limmad)[7]<-"entrez_id"
  colnames(limmad)[8]<-"Gene.symbol" #supporting Advaita
  colnames(limmad)[9]<-"ensembl_id"
  colnames(limmad)[10]<-"Gene.title" #supporting Advaita
  #grab the meta data matching the ensembl gene ids from limma
  Index <- mcols(rowRanges(kexp))$gene_id %in% limmad[,9]
  newFeatures <- mcols(rowRanges(kexp))[Index,]
  Features<-newFeatures[c(4,8:9)] #grabbing gene_id, gene_biotype and biotype_class
  uniqueFeatures<-Features[!duplicated(Features$gene_id),]
  limmad<-cbind(limmad,"NA")
  limmad<-cbind(limmad,"NA")
  colnames(limmad)[c(11:12)]<-c("gene_biotype","biotype_class")

  for(i in 1:nrow(limmad)) { # cbind biotype class to limma results
    indx <- which(rownames(limmad) == uniqueFeatures$gene_id[i])
    limmad[indx,c(11:12)] <- cbind(uniqueFeatures$gene_biotype[i],
                                   uniqueFeatures$biotype_class[i])
  }

  res$WithMeta<-limmad
  return(res)
} # }}} format limma results

