#' @title this analyzes downstream of wgcna
#' @description this can investigate annotations of the datTraits object plots individual color modules and an individual biotype correlation and pvalue significance 
#' @param datExpr list entry from wgcna part 1 method [[1]]
#' @param datTraits list entry from wgcna part 1 method [[[2]]
#' @param biotype   names of datTraits
wgcna_analysis<-function(datExpr,datTraits,biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"),biocolor="blue",whichWGCNA=c("single","block")){
  ##fix me: allow for multi color and biotype plots

load("wgcna.dataInput.RData")
###declare needed objects from load
message(paste0("found: ",names(lnames)))
bwModuleColors<-lnames["moduleColors"]
MEs<-lnames["MEs"]
datExpr<-lnames["datExpr"]
datTraits<-lnames["datTraits"]
annot<-lnames["annot"]

whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
biotype<-match.arg(biotype,c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"))
nGenes<-ncol(datExpr)
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
if(whichWGCNA=="single"){
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
} else {
MEs0 = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
}
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

####

weight<-as.data.frame(datTraits[,grep(biotype,colnames(datTraits))])
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)



geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

###
module<-biocolor
column<-match(module,modNames)
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",biotype),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

##map these to gene IDs
names(datExpr)[moduleColors==module]  
##need to form geneInfo0 using converted
  
  
stopifnot(nrow(geneTraitSignificance)==ncol(datExpr))
stopifnot(nrow(GSPvalue)==ncol(datExpr))
# Create the starting data frame


if(whichWGCNA=="single"){
geneInfo0 = data.frame(substanceBXH = names(datExpr),
                      geneSymbol = annot[probes2annot,]$mgi_symbol,
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)

} else {


geneInfo0 = data.frame(substanceBXH = names(datExpr),
                      geneSymbol = annot[probes2annot,]$mgi_symbol,
                      moduleColor = bwModuleColors,
                      geneTraitSignificance,
                      GSPvalue)


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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

return(geneInfo)
} #main
