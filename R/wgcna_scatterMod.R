#' @title plots a scatter correlation image of a fixed biotype and fixed eigengene module color
#' @description allows for exploration of data presented in wgcna_Cormap
#' @param lnames a list of the results from the call wgcna method
#' @param biocolor  the color of eigengene module
#' @param biotype   the tx biotype of interest
#' @import WGCNA
#' @export
#' @return images of plot
wgcna_scatterMod<-function(lnames,biocolor="blue",biotype="ERVL"){

if(is.null(lnames)==FALSE){
message("found wgcna data")
} else {
 stop("please run wgcna before analysis")
}
bwModuleColors<-lnames[["moduleColors"]]
MEs<-lnames[["MEs"]]
datExpr<-lnames[["datExpr"]]
datTraits<-lnames[["datTraits"]]
annot<-lnames[["annot"]]


#biotype<-match.arg(biotype,c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"))
nGenes<-ncol(datExpr)
nSamples = nrow(datExpr);
weight<-as.data.frame(datTraits[,grep(biotype,colnames(datTraits))])
names(weight) = as.character(biotype)

MEs0 = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# names (colors) of the modules
###
moduleColors<-bwModuleColors
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
readkey()


}
