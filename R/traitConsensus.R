#' @title finds consensus amongst desired traits
#' @description for traits with similiar correlation patters for a given module, this will find a consensus amongst groups
#' @import TxDbLite
#' @export
#' @return a df of the most significant common elements across biotypes per a given module
traitConsensus<-function(dbname=NULL,dbPath=".",consensus=c("Alu","ERVK","ERV3","ERVL","LTR Retrotransposon"),Module.color="blue",species=c("Homo.sapiens","Mus.musculus")){

  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  stopifnot(file.exists(paste0(dbPath,"/",dbname)))
  dbn<-paste0(dbPath,"/",dbname)

  colors<-dbListTables(dbconn(wgcnaDbLite(basename(dbn))))
  stopifnot(Module.color%in%colors)
  traitSet<-traitsBy(wgcnaDbLite(basename(dbn)),Module.color=Module.color)

  cons.id<-substring(names(traitSet),length(unlist(strsplit(Module.color,"",fixed=T)))+2)%in%consensus
   trait.cons<-traitSet[cons.id]
   intCons<-intersect(trait.cons[[1]]$row_names,trait.cons[[2]]$row_names)
    for(i in 3:length(names(trait.cons))){
    intCons<-intersect(intCons,trait.cons[[i]]$row_names)
    }

  ###FIX ME::: include the average of the p.values must report the p.values
  gxs<-data.frame(gene_id=intCons,stringsAsFactors=FALSE)  
  rownames(gxs)<-gxs$gene_id
   ##FIX ME: annotate intCons
    entrezID <- getEntrezIDs(gxs, species)
   entrez.match<-match(rownames(gxs),names(entrezID))
  stopifnot(all(rownames(gxs)==names(entrezID)[entrez.match]))
  gxs$entrezid<-entrezID[entrez.match]
  symbolNames<-getSymbols(gxs,species)
  symbol.match<-match(rownames(gxs),names(symbolNames))
  stopifnot(all(rownames(gxs)==names(symbolNames)[symbol.match]))
  gxs$hgnc_symbol<-symbolNames[symbol.match]
  write.csv(gxs,file=paste0(paste(consensus,collapse="."),"_",Module.color,".consensus.csv"),quote=FALSE,row.names=TRUE)
  return(gxs)

} ##main



#' @describeIn repeatTablesFromWGCNA
#' 
#' add EntrezGene IDs for Ensembl genes
#'
#' @param gxs       geneIDs from a dataframe
#' @param species  what kind of organism these genes are from
#' 
#' @return  entrez_id values for the genes, where found
#'
#'
getEntrezIDs <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)
  res <- try(mapIds(get(org$package), keys=as.character(gxs$gene_id),
                column="ENTREZID", keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No ENTREZID mappings were found for these genes...")
    return(rep(NA, length(gxs)))
  } else {
    return(res)
  }
} # }}}



#' @describeIn repeatTablesFromWGCNA
#' 
#' add symbols for Ensembl genes
#' @return  symbols for the genes, where found 
#'
#' 
getSymbols <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)

  ## needed since Ens83...
  if (any(grepl("\\.", gxs$gene_id))){
   gxs$gene_id<-gsub("\\.","",gxs$gene_id)
  }#no period in gene names.... 

  res <- try(mapIds(get(org$package), keys=as.character(gxs$gene_id),
                    column=org$symbol, keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No SYMBOLS were found for these genes...")
    return(rep(NA, length(gxs)))
  } else {
    return(res)
  }
} # }}}


