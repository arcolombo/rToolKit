#' @title GO enrichment for wgcna_analysis modules
#' @description downstream from wgcna_analysis we can investigate enrichment
wgcna_goEnrich<-function(lnames,intModules=c("blue","steelblue","darkturquoise"),species=c("human","mouse"),entrezOnly=TRUE ) {

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
intModules<-intColors
} else{
probes = names(datExpr)
id<-!is.na(annot$entrezgene)
probes2annot<-match(probes,annot$ensembl_gene_id)
allEntrezID<-annot$entrezgene[probes2annot]
intModules<-intColors
}


for (module in intModules)
{
  # Select module probes
  t<-showModuleMembers(lnames,biocolor=module)
 
  # Write them into a file
  fileName = paste("IDs-", module, ".txt", sep="");
  write.table(as.data.frame(t), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.

###FIX ME:  does bwModuleColors need to match allEntrez Id?
allEntrez<-annot$entrezgene
GOenr = GOenrichmentAnalysis(bwModuleColors, allEntrez, organism = species, nBestP = 15);

tab = GOenr$bestPTerms[[4]]$enrichment

write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

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
