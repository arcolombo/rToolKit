#' @title this will query specific genes in the filteredCorSet lists
#' @description this enriches each module to understand the function shown in Cormap creates a cross-corrleation data.frame and annotates it for each tile in teh Cormap. then uses the kexp to run a paired qusage enrichment analysis for each module. One should only enrich the genes in Module, irrespective of the tx.biotype, because the genes are fixed across all biotypes, so the enrichment patterns are identical.  so 1 and only 1 biotype needs to be interesting
#' @import qusage
#' @import WGCNA
#' @import arkas
#' @import edgeR

#' @return qusage data
wgcna_qusage<-function(kexp,lnames,geneModuleDF=NULL,biotype="Alu",intMods=c("blue","cyan","brown","purple","red","lightyellow","greenyellow","lightgreen","turquoise"),p.cutoff=0.07,MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),comparison="pHSC",control="LSC",how=c("cpm","tpm") ){
  ###THIS IS DEPRECATED most of this function is a weaker form of database creatoin now overpowered by wgcnaDbLite 
  how<-match.arg(how,c("cpm","tpm"))
   geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))
  if(how=="cpm"){
  ##TMM normalize
  counts<-collapseBundles(kexp,"gene_id") 
  dge<-DGEList(counts=counts)
  dge<-calcNormFactors(dge)
  expr<-cpm(dge,log=FALSE)
  counts<-log2(1+expr) ##enirhcment on log2 is required
  } else{
  counts<-collapseTpm(kexp,"gene_id")
  counts<-log2(1+counts)
  }
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
 
  weight<-data.frame(biotype=datTraits[,grep(biotype,colnames(datTraits))])
  colnames(weight)<-biotype
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

    qusage_run1<-qusageRun(cnts_mt=cnts_pL,MsigDB=MsigDB,comparison="pHSC",control="LSC",module=module)
    qusage_run2<-qusageRun(cnts_mt=cnts_bL,MsigDB=MsigDB,comparison="Blast",control="LSC",module=module)
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
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, ylim=c(-1,1))

   dev.off()
  } else{
   pdf(paste0("filtered_",p.cutoff,"_drivers_",module,".",colnames(geneTraitCor)[i],".pdf"))
   verboseScatterplot(df[,1],
                   df[, 2],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",colnames(geneTraitCor)[i]),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,ylim=c(-1,1))

    dev.off()

  }
   } #verbose all types
 }##across all colors


} ##main
