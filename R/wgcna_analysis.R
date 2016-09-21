#' @title this analyzes downstream of wgcna
#' @description this can investigate annotations of the datTraits object plots individual color modules and an individual biotype correlation and pvalue significance. the scatter plot returned is produced by the correlation value, p.value and regression line; it uses the base function lm to fit a linear model.
#' @param datExpr list entry from wgcna part 1 method [[1]]
#' @param datTraits list entry from wgcna part 1 method [[[2]]
#' @param biotype   names of datTraits
wgcna_analysis<-function(lnames,biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"),biocolor="blue",whichWGCNA=c("single","block"),intColors=c("blue","brown"),useBiCor=TRUE){
  ##fix me: allow for multi color and biotype plots
  ##FIX ME: load the MEs don't recalculate.  whichWGCNA shouldnt be used???
  if(is.null(lnames)==TRUE){
  load("wgcna.dataInput.RData")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  MEs<-lnames[["MEs"]]
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
  } else if(is.null(lnames)==FALSE){
  message("found wgcna.dataInput")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
   MEs<-lnames[["MEs"]]
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]

  } else {
  stop("please run wgcna, need the output before analyzing.")
  }


  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  biotype<-match.arg(biotype,c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"))
  nGenes<-ncol(datExpr)
  nSamples = nrow(datExpr);
 #plot adjacency heatmap of each biotype 
  wgcna_adjacencyHeatmap(MEs,datTraits)


# names (colors) of the modules
  modNames = substring(names(MEs), 3)
geneTraitModuleDF<-wgcna_scatterMod(lnames) ##plots all colors across all biotypes...  
##from examining scatter mod we have abs correlations with at least 0.50
## brown: Alu
## magenta: ERV1,ERVK,L1,End.Retr
## blue : ERVK, ERVL, L1
## salmon: ERV3, ERVK,ERVL, L1
## tan: End.Retro, ERV1,ERV3,ERVL
## grey60: End.Ret, ERVL
## black: Alu
## darkgrey: L2
## steelblue: ERVK, L1

##TO DO:
##screen and filter the geneTraitDF
##perhaps short list
##chooseTopHub
## fdr filter
##permutation test (screen)
## enrich using qusage
## call annotations to get top genes
## regression methods: cox regression, votingLinearPredictions  the verboseScatter returns Lm with really small pvalues...cant interpret.

#####
  probes<-names(datExpr)
  probes2annot<-match(probes,annot$ensembl_gene_id)
  if(whichWGCNA=="single"){
  geneInfo0 = data.frame(substanceBXH = names(datExpr),
                      geneSymbol = annot[probes2annot,]$mgi_symbol,
                      entrez = annot[probes2annot,]$entrezgene,
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)

} else {

##will have some NAs
geneInfo0 = data.frame(substanceBXH = names(datExpr),
                      geneSymbol = annot[probes2annot,]$mgi_symbol,
                      entrez = annot[probes2annot,]$entrezgene,
                      moduleColor = bwModuleColors[probes2annot],
                      geneTraitSignificance[probes2annot,],
                      GSPvalue[probes2annot,])


} #block level


# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$geneTraitSignificance.probes2annot...));
geneInfo = geneInfo0[geneOrder, ]

######plots final ##################################
probes = names(datExpr)
probes2annot<-match(probes,annot$ensembl_gene_id)
allEntrezID<-annot$entrezgene[probes2annot] 
intModules<-intColors

for (module in intModules)
{
  # Select module probes
  t<-showModuleMembers(lnames,biocolor=module)
 
  # Write them into a file
  fileName = paste("EntrezIDs-", module, ".txt", sep="");
  write.table(as.data.frame(t), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.

############GO
##FIX ME:::: add qusage here , this method is not recommended
allEntrezID<-allEntrezID[!is.na(allEntrezID)] ###FILTER OUT NA ENTREZ

GOenr = GOenrichmentAnalysis(moduleColors, allEntrezID, organism = "mouse", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
###QUSAGE HERE

################

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
screenTab



return(geneInfo)
} #main
