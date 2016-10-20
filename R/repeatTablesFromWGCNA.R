#' @title creates a WGCNA annotation package that includes many tables
#' @description  after WGCNA creates a data object post networking and partitioning of data tables can add structure to the many gene~module and gene~trait correlations both including student-t and exact Fisher p.values with Gene annotations.  The tables created will be one table for each gene~module collection, and one table for each gene~trait collection with annotations. NOTE: wgcna.R needs to be called and the network object needs to be created prior the database tables.Further if you are not building a data base with phenotypic traits, then the datTraits will be NA.
#' @param kexp any kexp doesn't matter, only use the metadata
#' @param lnames is the saved object from wgcna.R that includes basic formatted data compatible with WGCNA, datExpr, datTraits, annot library etc. Note that datTraits will be NA if considering gene modules only.
#' @param useBiCor boolean ,  if TRUE uses mid-correlation mid-weight robust to outliers
#' @param verbose boolean true for std out
#' @param dbname  character
#' @param annotate boolean, for genes level set to TRUE, repeats set to FALSE
#' @param version character for meta
#' @param byWhich character 
#' @importFrom Rsamtools indexFa index
#' @importFrom Rsamtools scanFaIndex FaFile scanFa path
#' @importFrom GenomeInfoDb seqnames
#' @importFrom DBI dbConnect dbDriver dbWriteTable dbGetQuery dbDisconnect
#' @importFrom S4Vectors DataFrame
#' @export
#' @return a SQLlite Data base 
repeatTablesFromWGCNA<-function(kexp,lnames,useBiCor=TRUE,verbose=TRUE,dbname="wgcnaDBLite",annotate=FALSE,version="1.0.0",byWhich=c("gene","repeat") ){
  
   byWhich<-match.arg(byWhich,c("gene","repeat"))
   tnxomes<-transcriptomes(kexp)
   ##supporting only ENSEMBL
   tnxomes<-strsplit(tnxomes,",")
   Ensmbl<-unlist(tnxomes)[grep("Ens",unlist(tnxomes))]
  
 if(verbose) cat("Extracting Networking data...\n")
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  sat.ID<-grep("satellite",colnames(datTraits),ignore.case=T)
  if(length(sat.ID)>1){ 
  colnames(datTraits)[sat.ID[1]]<-"satelliteBB"
  }
  how<-lnames[["how"]]
  if(annotate==FALSE){
  #if annotate is FALSE, then we use the default annotation data object, however the colnames need to be changed from lnames ensembl_gene_id to gene_id (crucual for uniformaity with dbLite.
  annot<-lnames[["annot"]]
  yy<-grep("ensembl_gene_id",colnames(annot))
  colnames(annot)[yy]<-"gene_id"
  }
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
  modNames = substring(names(MEs), 3)
  nGenes<-ncol(datExpr)
  nSamples = nrow(datExpr);
  modNames<-substring(colnames(MEs),3)
  moduleColors<-bwModuleColors

  
  weight<-as.data.frame(datTraits)
  
  if(byWhich=="gene"){ 
  if(any(grepl("^ENSG",colnames(datExpr)))){
   species<-"Homo.sapiens"
   packageName<-gsub("EnsDbLite","wgcnaDbLite",Ensmbl)
  } else if(any(grepl("^ENSMUSG",colnames(datExpr)))){
  species<-"Mus.musculus"
  packageName<-gsub("EnsDbLite","wgcanDbLite",Ensmbl)
  }
 } else{
 cat("Homo.Sapiens supported for repeat tables...\n")
    species<-"Homo.sapiens"
   packageName<-gsub("EnsDbLite","wgcnaDbLite",Ensmbl)
 }

 if(useBiCor==FALSE){
  ##finds correlation per gene in modules
 if(verbose) cat("Calculating Module Correlations with EigenGenes\n")
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ##pvalue per gene per module
  if(verbose) cat("Calculating student-t pvalue ...\n")
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p_MM", modNames, sep="");
  if(verbose) cat("Calculating Cross Correlation of genes~traits...\n")
  geneTraitCor = as.data.frame(cor(datExpr, weight, use = "p"));
  if(verbose) cat("Calculating student-t pvalues of gene~trait correlations...\n")
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
  if(verbose) cat("Calculating Fisher-Exact p.values...\n")  
  geneTraitFisherPvalue = as.data.frame(corPvalueFisher(as.matrix(geneTraitCor), nSamples));
  colnames(geneTraitCor)<-paste("GCor_",colnames(weight),sep="");
  colnames(geneTraitPvalue)<-paste("p_GCor_",colnames(weight),sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf_GCr_",colnames(weight),sep="");
  } else {
  if(verbose) cat("Calculating Module bi-Correlations with EigenGenes...\n")
 geneModuleMembership = as.data.frame(bicor(datExpr, MEs, use = "all.obs"));
  if(verbose) cat("Calculating student-t p.values...\n")
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); ##pvalue per gene in each module
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p_MM", modNames, sep="");
  ##gene pvale per trait
  if(verbose) cat("Calculating bi-correlation of gene~trait...\n")
  geneTraitCor = as.data.frame(bicor(datExpr, weight, use = "all.obs"));
  if(verbose) cat("Calculating student-t p.values gene~trait...\n")
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
  geneTraitFisherPvalue<-as.data.frame(corPvalueFisher(as.matrix(geneTraitCor),nSamples));
   colnames(geneTraitCor) = paste("GCor_", colnames(weight), sep="");
  colnames(geneTraitPvalue) = paste("p_GCor_", colnames(weight), sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf_GCr_",colnames(weight),sep="");
  }

  stopifnot(all(rownames(geneModuleMembership)==rownames(MMPvalue))==TRUE)
  geneModuleDF<-cbind(geneModuleMembership,MMPvalue)
  stopifnot(all(rownames(geneModuleDF)==rownames(geneTraitCor))==TRUE)
  geneModuleDF<-cbind(geneModuleDF,geneTraitCor)
  stopifnot(all(rownames(geneModuleDF)==rownames(geneTraitPvalue))==TRUE)
  geneModuleDF<-cbind(geneModuleDF,geneTraitPvalue)
  geneModuleDF<-cbind(geneModuleDF,geneTraitFisherPvalue)
 
  ##call annotations?
   if(annotate==TRUE){
  annot<-data.frame(gene_id=rownames(geneModuleDF),row.names=rownames(geneModuleDF))
  entrezID <- getEntrezIDs(annot, species)
  entrez.match<-match(rownames(annot),names(entrezID))
  stopifnot(all(rownames(annot)==names(entrezID)[entrez.match]))
  annot$entrezid<-entrezID[entrez.match]
  symbolNames<-getSymbols(annot,species)
  symbol.match<-match(rownames(annot),names(symbolNames))
  stopifnot(all(rownames(annot)==names(symbolNames)[symbol.match]))
  annot$hgnc_symbol<-symbolNames[symbol.match]
  geneModuleDF<-cbind(geneModuleDF,annot)
 } else{
  ##do not annotate namely for Repeat Elements
  if(byWhich=="repeat"){
  #for repeat we want the description in the annotation to be the hgnc_symbol, the description for wrcna is the tx biotype. more useful here
  annotDF<-annot[rownames(annot)%in%rownames(geneModuleDF),]
  gene.ID<-grep("gene_id",colnames(annotDF))
  ent.ID<-grep("entrezgene",colnames(annotDF))
  desc.ID<-grep("description",colnames(annotDF))
  annotDF<-annotDF[,c(gene.ID,ent.ID,desc.ID)]
  colnames(annotDF)<-c("gene_id","entrezid","hgnc_symbol")
  entrez.match<-match(rownames(annotDF),rownames(geneModuleDF))
  stopifnot(all(rownames(annotDF)==rownames(geneModuleDF[entrez.match,])))
  geneModuleDF<-cbind(geneModuleDF,annotDF)
  } else { 
  stop("For gene repeat Tables please flag annotate as TRUE...\n")
  }
 } ##annotate as FALSE



 ##the goal is to have a SQL library instead of a large data frame geneTraitModuleDF
 ##the bwLabels is the key, and should be written to a table as well. this is the cornerstone to accessing moduleGenes and all associated data.

 if(verbose) cat("Creating Database...\n")
  
  dbname<-paste(dbname,"sqlite",sep=".")
  con<-dbConnect(dbDriver("SQLite"),dbname=dbname)
  if(verbose) cat("Extracting Module Color Names...")
   modulesNeeded<-unique(moduleColors)
   modulePvalue.Needed<-colnames(geneModuleDF)[grep("p_MM",colnames(geneModuleDF))]
  stopifnot(length(modulesNeeded)==length(modulePvalue.Needed))
  if(verbose) cat("Extracting Pheno-types of interest...\n")
  traitsNeeded<-colnames(datTraits)
  geneTraitCorNeeded<-paste0("GCor_",traitsNeeded)
  pTraitNeeded<-paste0("p_GCor_",traitsNeeded)
  pfTraitNeeded<-paste0("pf_GCr_",traitsNeeded) 
  stopifnot(all(geneTraitCorNeeded%in%colnames(geneModuleDF)))
  stopifnot(all(pTraitNeeded %in% colnames(geneModuleDF)))
  stopifnot(all(pfTraitNeeded%in%colnames(geneModuleDF)))
  
  geneTraitCor.id<-match(geneTraitCorNeeded,colnames(geneModuleDF))
  pTrait.id<-match(pTraitNeeded,colnames(geneModuleDF))
  pfTrait.id<-match(pfTraitNeeded,colnames(geneModuleDF))
  #color,GCor_TxBio,pvalue , pf.value,annotations
  ## FIX ME: add annotations
  annot.id<-match(rownames(geneModuleDF),annot$gene_id)
  
  for(i in 1:length(modulesNeeded)){
  mod.id<-which(modulesNeeded[i]==substring(colnames(geneModuleDF),3 ))
  mod.p.id<-which(modulePvalue.Needed[i]==substring(colnames(geneModuleDF),5))
  colorTable<-geneModuleDF[,c(mod.id,mod.p.id,geneTraitCor.id,pTrait.id,pfTrait.id)]
  colorKey<-as.factor(moduleColors)
  colorTable<-cbind(colorTable,colorKey)
  colorTable<-cbind(colorTable,annot[annot.id,]) 
  stopifnot(all(rownames(colorTable)==rownames(geneModuleDF)))
  stopifnot(all(rownames(colorTable)==colorTable$gene_id))
  
  ##write dbLite table  
  if(verbose) cat(paste0("Creating the database for ",modulesNeeded[i],"\n"))
  dbWriteTable(con,name=modulesNeeded[i],colorTable,overwrite=T,row.names=T)
  dbGetQuery(con,paste0("create index colorKey_idx","_",i," on ",modulesNeeded[i]," (colorKey);"))
   }
  ##FIX ME :  write out some meta data color module summaries, species .. 
  Metadata <- wgcnaDbLiteMetadata(kexp,
                               packageName=packageName,
                               species=species,
                               Ensmbl=Ensmbl,
                               how=how,
                               byWhich=byWhich)
 dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)




  cat("done. \n")
 dbDisconnect(con)
 return(dbname)
} #main



#' @describeIn repeatTablesFromWGCNA
#' 
#' add EntrezGene IDs for Ensembl genes
#'
#' @param gxs       geneIDs from a dataframe
#' @param species  what kind of organism these genes are from
#' @export 
#' @return  entrez_id values for the genes, where found
#'
#'
getEntrezIDs <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)
  res <- try(mapIds(get(org$package), keys=as.character(gxs$gene_id),
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
#' @export
#' @return  symbols for the genes, where found 
#'
#' 
getSymbols <- function(gxs, species) { # {{{
  org <- getOrgDetails(species)
  library(org$package, character.only=TRUE)

  ## needed since Ens83...
  if (any(grepl("\\.", gxs$gene_id))){
   gxs$gene_id<-gsub("\\.","",gxs$gene_id)
  }#no period in gene names.... 

  res <- try(mapIds(get(org$package), keys=as.character(gxs$gene_id),
                    column=org$symbol, keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No SYMBOLS were found for these genes...")
    return(rep(NA, length(gxs)))
  } else {
    return(res)
  }
} # }}}

#' @describeIn repeatTablesFromWGCNA 
#' 
#' create metadata for an wgcnaDbLite instance
#'
#' @param packageName   the name of the annotation package to be built 
#' @param genomeVersion name of genome assembly for coordinates, e.g. "GRCh38"
#' @param sourceFile    name of FASTA file(s) whence it was built, as a string
#' 
#' @return a data.frame of metadata suitable for cramming into the database
#'
#' @export
wgcnaDbLiteMetadata <- function(kexp,packageName,species=NULL,Ensmbl,how=how,byWhich=byWhich) { # {{{

  
  Ensmbl<-gsub(" ","",Ensmbl)
  org <- getOrgDetails(species)
  sourceFile <- Ensmbl
  kversion<-kallistoVersion(kexp)
  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "wgcnaDbLite")
  MetaData[3,] <- c("type_of_gene_id", "Ensembl Gene ID")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", species)
  MetaData[7,] <- c("kallistoVersion", kversion)
  MetaData[8,] <- c("sourceFile", sourceFile)
  MetaData[9,]<- c("how",how)
  MetaData[10,]<-c("which",byWhich)
  return(MetaData)

} # }}}

