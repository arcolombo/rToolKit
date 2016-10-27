#' @title GO enrichment for wgcna_analysis modules
#' @description downstream from wgcna_analysis we can investigate enrichment. The wgcna_goEnrich is used when creating the repeatTables SQL database.  the WGCNA:goEnrich is a bit slow, so it is ideal to call this method once upon SQLdb creation and query the info needed.  so repeatTablesFromWGCNA calls this method and stores into the DB.  WGCNA:goEnrich is an experimental method, and repeatToolKit supports qusage, however qusage is differential activity, where goEnrich tests for overall activity as a preliminary glimpse.  useful for exploration
#' @param lnames  this is the main method from wgcna, this is a slow calculation which sucks, but you call it once, and save it, so you can do downstream analysis with one painful object creation call
#' @param species human or mouse for annotation
#' @param entrezOnly boolean, I am not sure if you should use all entries for goEnrichment, with and without entrez.  so keeping a flag for now.  it seems as though the GoEnirch call filters NA entrez automatically. so for safe keeping keep entrezOnly to false, and let WGCNA enrichment caller filter out the NAs.
#' @export
#' @return images and go enrichment object
wgcna_goEnrich<-function(lnames,species=c("human","mouse"),entrezOnly=FALSE ) {

  species<-match.arg(species,c("human","mouse"))
  if(is.null(lnames)==TRUE){
  load("wgcna.dataInput.RData")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  } else if(is.null(lnames)==FALSE){
  message("found wgcna.dataInput")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  } else {
  stop("please run wgcna, need the output before analyzing.")
  }


###need MEs, moduleLables,moduleColors, and geneTree loaded
  if(entrezOnly==FALSE){
  probes = names(datExpr)
  probes2annot<-match(probes,annot$ensembl_gene_id)
  allEntrezID<-annot$entrezgene[probes2annot]
 } else{
  probes = names(datExpr)
  id<-!is.na(annot$entrezgene)
  annot2<-annot[id,] ##short list to valid entrezID
  probes2annot<-match(probes,annot2$ensembl_gene_id) #no NA entrez matches
  stopifnot(annot2$esembl_gene_id[probes2annot]==probes)
  allEntrezID<-annot2$entrezgene[probes2annot]
  }


# As background in the enrichment analysis, we will use all probes in the analysis.
###FIX ME:  examine goEnrichment script
  GOenr = GOenrichmentAnalysis(bwModuleColors, 
                               allEntrezID, 
                               organism = species, 
                               nBestP = 15);

  tab = GOenr$bestPTerms[[4]]$enrichment
  write.table(tab,
              file = "GOEnrichmentTable.csv", 
              sep = ",", 
              quote = TRUE, 
              row.names = FALSE)

  keepCols = c(1, 2, 5, 6, 7, 12, 13);
  screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
  numCols = c(3, 4);
  screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
  screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
  colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
  rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
  options(width=95)
# Finally, display the enrichment table:
  print(screenTab)
  return(screenTab)

} #main
