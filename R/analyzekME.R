#' @title analyzes kMEs in comparison and control MEs
#' @description given a control module color set we can test for conservation in a comparison set.  by calculating the signedkME, which is the connectivity of the ModuleEigen values of the first PC1, this will quantify the correlation of the gene and moduleEigenvalues between modules across data sets.
analyzekME<-function(datExpr1=NULL,datExpr2=NULL,ME_G1=NULL,colorsG1=NULL,modulesG1=NULL,species="Homo.sapiens",main="Gene"){
  main<-match.arg(main,c("Gene","Repeat"))
 geneModuleMembership1 = signedKME(t(datExpr1), ME_G1)
  colnames(geneModuleMembership1)=paste("PC",colorsG1,".cor",sep="");
  MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExpr1)[[2]]);
  colnames(MMPvalue1)=paste("PC",colorsG1,".pval",sep="");
  Gene= rownames(datExpr1)
  kMEtable1 = cbind(Gene,Gene,modulesG1)
  for (i in 1:length(colorsG1)){
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
  }
 colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1),
  colnames(MMPvalue1))))

  write.csv(kMEtable1,file=paste0(main,"_kMEtable1.csv"),row.names=FALSE)


 # First calculate MEs for A2, since we haven't done that yet with resepect to Modules in  control
  PCs2g = moduleEigengenes(t(datExpr2), colors=modulesG1)
  ME_2g = PCs2g$eigengenes
  geneModuleMembership2 = signedKME(t(datExpr2), ME_2g)
  colnames(geneModuleMembership2)=paste("PC",colorsG1,".cor",sep="");
  MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(datExpr2)[[2]]);
  colnames(MMPvalue2)=paste("PC",colorsG1,".pval",sep="");
  kMEtable2 = cbind(Gene,Gene,modulesG1)
  for (i in 1:length(colorsG1)){
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) }
  colnames(kMEtable2)=colnames(kMEtable1)

  write.csv(kMEtable2,file=paste0(main,"_kMEtable2.csv"),row.names=FALSE)

par(mfrow=c(1,2),mar=c(5,4,4,2))
for (c in 1:length(colorsG1)){
verboseScatterplot(geneModuleMembership2[,c],
                   geneModuleMembership1[,c],
                   main=colorsG1[c],
                   cex.axis=0.8,
                   cex.main=0.9,
                   cex.lab=1,
                  xlab="kME in A2",ylab="kME in A1")
readkey()
}; 

   for (c in 1:length(colorsG1)){
   inMod = modulesG1== colorsG1[c]
   verboseScatterplot(geneModuleMembership2[inMod,c],
                      geneModuleMembership1[inMod,c],
                      main=colorsG1[c],
                      xlab="kME in A2",ylab="kME in A1")
     readkey()   
   };
  ##produce images.
  pdf(paste0(main,"_all_kMEtable2_vs_kMEtable1.pdf"),height=8,width=8)
  par(mfrow=c(1,2),mar=c(5,4,4,2))
  for (c in 1:length(colorsG1)){
  verboseScatterplot(geneModuleMembership2[,c],
                   geneModuleMembership1[,c],
                   main=colorsA1[c],
                   xlab="kME in A2",
                   ylab="kME in A1")
  };
   dev.off()
  pdf(paste0(main,"_inModule_kMEtable2_vs_kMEtable1.pdf"),height=8,width=8)
  par(mfrow=c(1,2),mar=c(5,4,4,2))
  for (c in 1:length(colorsG1)){
  inMod = modulesG1== colorsG1[c]
  verboseScatterplot(geneModuleMembership2[inMod,c],
                     geneModuleMembership1[inMod,c],
                     main=colorsA1[c],
                     xlab="kME in A2",ylab="kME in A1")
   }; 
  dev.off()


  topGenesKME = NULL
  for (c in 1:length(colorsG1)){
  kMErank1= rank(-geneModuleMembership1[,c])
  kMErank2= rank(-geneModuleMembership2[,c])
  maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=25])
  }; 
  colnames(topGenesKME) = colorsG1
   topGenesKME<-as.data.frame(topGenesKME)
    symbolsKME<-NULL
  ##annotate!!!
   if(main=="Gene"){
    for(i in 1:ncol(topGenesKME)){
   symbolNames<-getSymbols(topGenesKME[,i],species)
   symbolsKME<-cbind(symbolsKME,symbolNames)
  }
  colnames(symbolsKME)<-colnames(topGenesKME)
   } else{
   symbolsKME<-topGenesKME
   }
    cat("hubs, which are significantly correlated to ME and high within module connectivity, found in both gene networks with the highest kME values:\n")
   print(symbolsKME)
return(symbolsKME)

} ##main



#' @describeIn repeatTablesFromWGCNA
#' 
#' add EntrezGene IDs for Ensembl genes
#'
#' @param gxs       geneIDs from a dataframe
#' @param species  what kind of organism these genes are from
#' 
#' @return  entrez_id values for the genes, where found
#'
#' 
getEntrezIDs <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)
  res <- try(mapIds(get(org$package), keys=as.character(gxs),
                column="ENTREZID", keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No ENTREZID mappings were found for these genes...")
    return(rep(NA, length(gxs)))
  } else {
    return(res)
  }
} # }}}




#' @describeIn repeatTablesFromWGCNA
#' 
#' add symbols for Ensembl genes
#' @return  symbols for the genes, where found 
#'
#'
getSymbols <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)

  ## needed since Ens83...
  if (any(grepl("\\.", gxs))){
   gxs<-gsub("\\.","",gxs)
  }#no period in gene names.... 

  res <- try(mapIds(get(org$package), keys=as.character(gxs),
                    column=org$symbol, keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No SYMBOLS were found for these genes...")
    return(rep(NA, length(gxs)))
  } else {
    return(res)
  }
} # }}}

