#' @title GO enrichment for wgcna_analysis modules
#' @description downstream from wgcna_analysis we can investigate enrichment
wgcna_goEnrich<-function(datExpr,datTraits,annot,intModules=c("blue","steelblue","darkturquoise") ) {

###need MEs, moduleLables,moduleColors, and geneTree loaded
probes<-names(datExpr)
probes2annot<-match(probes,annot$ensembl_gene_id)
allEnsIDs = annot$mgi_symbol[probes2annot];
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modEnsIDs = allEnsIDs[modGenes];
  # Write them into a file
  fileName = paste("IDs-", module, ".txt", sep="");
  write.table(as.data.frame(modEnsIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("EnsGIDs-all.txt", sep="");
write.table(as.data.frame(allEnsIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

###FIX ME:  GO requires entrez ID
GOenr = GOenrichmentAnalysis(bwModuleColors, names(datExpr), organism = "mouse", nBestP = 10);

} #main
