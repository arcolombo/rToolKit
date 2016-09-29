#' @title this will query specific genes in the filteredCorSet lists
#' @description this enriches each module to understand the function shown in Cormap creates a cross-corrleation data.frame and annotates it for each tile in teh Cormap. then uses the kexp to run a paired qusage enrichment analysis for each module tile.
#' @import qusage
#' @import WGCNA
#' @import arkas
#' @export
#' @return qusage data
wgcna_qusage<-function(kexp,lnames,geneModuleDF=NULL,biotype=c("ERV1","ERV2","Endogenous Retrovirus", "ERV3","ERVL", "L1","L2","LTR Retrotransposon"),intMods=c("blue","cyan","brown","purple","tan","red","lightyellow","greenyellow","lightgreen","darkturquoise","turquoise"),p.cutoff=0.07,MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),comparison="pHSC",control="LSC" ){

  ##The purpose of this function is to set up a control in a sense
  ## we have all of the modules and patterns across biotypes
  ## so we look at NK and ISGs and create a subset of the modules for known repeat associations.  then we can look at enrichment and DE on the entire module, or query DEs in each module.  we want to find out which module does the ISGs fall into ? which module does NK fall into?  DE ? 
   geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))
  counts<-collapseBundles(amlX,"gene_id")
  counts<-log2(1+counts) ##enirhcment on log2 is required
  if(is.null(geneModuleDF)==TRUE){
  geneModuleDF<-wgcna_scatterMod(lnames)
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
  modNames<-substring(colnames(MEs),3)
 moduleColors<-bwModuleColors
 intMods<-c("blue","cyan","brown","purple","tan","red","lightyellow",
   "greenyellow","lightgreen","darkturquoise","turquoise")

  weight<-as.data.frame(datTraits[,biotype%in%colnames(datTraits)])
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
  stopifnot(all(colnames(geneTraitCor)==substring(colnames(geneTraitPvalue),3)))

  for(j in 1:length(as.vector(unique(intMods)))){
  module<-as.vector(unique(intMods))[j]
  column<-match(module,modNames)
  moduleGenes = moduleColors==module;
  for(i in 1:ncol(geneTraitCor)){
  df<-as.data.frame((geneModuleMembership[moduleGenes, column]),
                   row.names=rownames(geneModuleMembership)[which(moduleGenes=="TRUE")])
  df<-cbind(df,geneTraitCor[moduleGenes,i])
  pvalue<-geneTraitPvalue[moduleGenes,i]
  df<-cbind(df,pvalue)
  df<-df[order(df$pvalue,decreasing=FALSE),]
  colnames(df)<-c(paste0("MM.",as.vector(unique(intMods))[j]),colnames(geneTraitCor)[i],"p.value")
  ##now annotate
  id<-match(rownames(df),annot$ensembl_gene_id)
  df<-cbind(df,annot[id,])
  ##run qusage 
  message(colnames(geneTraitCor)[i])
  ###qusage on non-filter module
  counts_id<-match(rownames(df),rownames(counts))
  counts_filter<-counts[counts_id,]
  stopifnot(all(rownames(df)==rownames(counts_filter)))
  ##must make two group contrasts pHSC-LSC, Blast-LSC
  pHSC.id<-grep("pHSC_",colnames(counts_filter))
  LSC.id<-grep("LSC_",colnames(counts_filter))
  Blast.id<-grep("Blast_",colnames(counts_filter))
  cnts_pL<-counts_filter[,c(pHSC.id,LSC.id)]
  cnts_bL<-counts_filter[,c(Blast.id,LSC.id)]
  annot.pL<-match(rownames(cnts_pL),annot$ensembl_gene_id)
  stopifnot(all(rownames(cnts_pL)==annot$ensembl_gene_id[annot.pL]))
   rownames(cnts_pL)<-annot$hgnc_symbol[annot.pL]
   ##puts hgnc as rownames, there will be some duplicates because ENSgIDs corrrespond to multiples
   annot.bL<-match(rownames(cnts_bL),annot$ensembl_gene_id)
   stopifnot(all(rownames(cnts_bL)==annot$ensembl_gene_id[annot.bL]))
   rownames(cnts_bL)<-annot$hgnc_symbol[annot.bL]

    qusage_run1<-qusageRun(cnts_mt=cnts_pL,MsigDB=MsigDB,comparison="pHSC",control="LSC",module=module,tx.biotype=colnames(geneTraitCor)[i])
    qusage_run2<-qusageRun(cnts_mt=cnts_bL,MsigDB=MsigDB,comparison="Blast",control="LSC",module=module,tx.biotype=colnames(geneTraitCor)[i])
    write.csv(qusage_run1,file=paste0("pHSC.v.LSC","_",module,".",colnames(geneTraitCor)[i],"_",MsigDB,".csv"),quote=FALSE,row.names=TRUE)
    write.csv(qusage_run2,file=paste0("Blast.v.LSC","_",module,".",colnames(geneTraitCor)[i],"_",MsigDB,".csv"),quote=FALSE,row.names=TRUE)

   ###save the module cross correlation with biotype object
   save(df,file=paste0(intMods[j],"_",colnames(geneTraitCor)[i],".RData"),compress=TRUE)
    df.filtered<-df[which(df$p.value<=0.07),] ##drivers!!! most exterme values in correlations most significantly correlated items.
   
   write.csv(df.filtered,file=paste0(intMods[j],"_drivers_",colnames(geneTraitCor)[i],".csv"))

  if(nrow(df.filtered)>10){
  pdf(paste0("filtered_",p.cutoff,"_drivers_",module,".",colnames(geneTraitCor)[i],".pdf"))
   verboseScatterplot(df[,1],
                   df[, 2],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",colnames(geneTraitCor)[i]),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
verboseScatterplot(df[which(df$p.value<=p.cutoff),1],
                   df[which(df$p.value<=p.cutoff), 2],
                   xlab = paste("Module Drivers in", module, "module"),
                   ylab = paste0("Gene significance ",p.cutoff," ",colnames(geneTraitCor)[i]),
                   main = paste("Module Drivers and Biotype Drivers\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

   dev.off()
  } else{
   pdf(paste0("filtered_",p.cutoff,"_drivers_",module,".",colnames(geneTraitCor)[i],".pdf"))
   verboseScatterplot(df[,1],
                   df[, 2],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",colnames(geneTraitCor)[i]),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

    dev.off()

  }
   } #verbose all types
 }##across all colors


} ##main
