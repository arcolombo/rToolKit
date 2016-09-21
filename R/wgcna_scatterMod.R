#' @title plots a scatter correlation image of a fixed biotype and fixed eigengene module color
#' @description allows for exploration of data presented in wgcna_Cormap
#' @param lnames a list of the results from the call wgcna method
#' @param biocolor  the color of eigengene module
#' @param biotype   the tx biotype of interest
#' @import WGCNA
#' @export
#' @return images of plot
wgcna_scatterMod<-function(lnames,biocolor="blue",biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"), useBiCor=TRUE){
  #FIX ME: the scatter correlation is not matching correlation from Cormap

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
  MEs<-lnames[["MEs"]]
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]

 nGenes<-ncol(datExpr)
  nSamples = nrow(datExpr);

###using all biotypes
 weight<-as.data.frame(datTraits[,biotype%in%colnames(datTraits)]

 if(useBiCor==FALSE){
  ##finds correlation per gene in modules
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ##pvalue per gene per module
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitCor = as.data.frame(cor(datExpr, weight, use = "p"));
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
 #  names(geneTraitPvalue) = paste("GS.", names(weight), sep="");
 # names(GSPvalue) = paste("p.GS.", names(weight), sep="");
} else{

 geneModuleMembership = as.data.frame(bicor(datExpr, MEs, use = "all.obs"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); ##pvalue per gene in each module
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  ##gene pvale per trait
  geneTraitCor = as.data.frame(bicor(datExpr, weight, use = "all.obs"));
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  # names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
  #names(GSPvalue) = paste("p.GS.", names(weight), sep="");
}
###FIX ME::: loop through biotypes per fix color.  hell add a double for loop
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
stopifnot(nrow(geneTraitSignificance)==ncol(datExpr))
stopifnot(nrow(GSPvalue)==ncol(datExpr))





######


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
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
readkey()


}
