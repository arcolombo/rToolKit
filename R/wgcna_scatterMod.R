#' @title plots a scatter correlation image of a fixed biotype and fixed eigengene module color
#' @description allows for exploration of data presented in wgcna_Cormap of individual modules plotted by individual traits and shows their absolute correlation.  also returns a data frame of modules' correlation and student and fisher pvalues in terms of ENSG id.  wgcna_filter is downstream of this analysis. note: you only want to call this function with plotALL=TRUE once, this will identify modules that have high absolute correlation. 
#' @param lnames a list of the results from the call wgcna method
#' @param biotype   the tx biotype of interest
#' @param useBiCor bicor is a WGCNA functoin that uses biweight midcorrelations and is robust against outliers
#' @import WGCNA
#' @export
#' @return prints images and returns module significance using student t test and fisher exact test
wgcna_scatterMod<-function(lnames,biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"), useBiCor=TRUE){
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
  modNames = substring(names(MEs), 3)
  nGenes<-ncol(datExpr)
  nSamples = nrow(datExpr);

###using all biotypes
 weight<-as.data.frame(datTraits[,biotype%in%colnames(datTraits)])

 if(useBiCor==FALSE){
  ##finds correlation per gene in modules
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ##pvalue per gene per module
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitCor = as.data.frame(cor(datExpr, weight, use = "p"));
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
   geneTraitFisherPvalue = as.data.frame(corPvalueFisher(as.matrix(geneTraitCor), nSamples));
  colnames(geneTraitCor)<-paste("GCor.",colnames(weight),sep="");
  colnames(geneTraitPvalue)<-paste("p.GCor.",colnames(weight),sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf.GCr.",colnames(weight),sep="");
  } else {

 geneModuleMembership = as.data.frame(bicor(datExpr, MEs, use = "all.obs"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); ##pvalue per gene in each module
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  ##gene pvale per trait
  geneTraitCor = as.data.frame(bicor(datExpr, weight, use = "all.obs"));
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
  geneTraitFisherPvalue<-as.data.frame(corPvalueFisher(as.matrix(geneTraitCor),nSamples));
   colnames(geneTraitCor) = paste("GCor.", colnames(weight), sep="");
  colnames(geneTraitPvalue) = paste("p.GCor.", colnames(weight), sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf.GCr.",colnames(weight),sep="");
  }

  cat("Plotting module and gene correlations...\n")
  moduleColors<-bwModuleColors
  pdf("scatterModulePlots.pdf")
  for(j in 1:length(as.vector(unique(moduleColors)))){

  module<-as.vector(unique(moduleColors))[j]
  column<-match(module,modNames)
  moduleGenes = moduleColors==module;

  for(i in 1:ncol(geneTraitCor)){
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot((geneModuleMembership[moduleGenes, column]),
                   (geneTraitCor[moduleGenes, i]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",colnames(geneTraitCor)[i]),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  #readkey()
   } #verbose all types
 }##across all colors

######FIX ME add a MSE plot
 for(j in 1:length(as.vector(unique(moduleColors)))){

  module<-as.vector(unique(moduleColors))[j]
  column<-match(module,modNames)
  moduleGenes = moduleColors==module;

  for(i in 1:ncol(geneTraitCor)){
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseIplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitCor[moduleGenes, i]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",colnames(geneTraitCor)[i]),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
 # readkey()
   } #verbose all types
 }##across all colors
dev.off()

####



  stopifnot(all(rownames(geneModuleMembership)==rownames(MMPvalue))==TRUE)
  geneModuleDF<-cbind(geneModuleMembership,MMPvalue)
  stopifnot(all(rownames(geneModuleDF)==rownames(geneTraitCor))==TRUE)
  geneModuleDF<-cbind(geneModuleDF,geneTraitCor)
  stopifnot(all(rownames(geneModuleDF)==rownames(geneTraitPvalue))==TRUE)
  geneModuleDF<-cbind(geneModuleDF,geneTraitPvalue)
  geneModuleDF<-cbind(geneModuleDF,geneTraitFisherPvalue)
  return(geneModuleDF)
}
