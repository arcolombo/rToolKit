#' @title creates a WGCNA annotation package that includes many tables
#' @description  after WGCNA creates a data object post networking and partitioning of data tables can add structure to the many gene~module and gene~trait correlations both including student-t and exact Fisher p.values with Gene annotations.  The tables created will be one table for each gene~module collection, and one table for each gene~trait collection with annotations.   FIX ME: there must be accessing methods for these tables.
#' @param lnames is the saved object from wgcna.R that includes basic formatted data compatible with WGCNA, datExpr, datTraits, annot library etc. FIX ME: have this annotation library be independent of wgcna.R
#' @param useBiCor boolean ,  if TRUE uses mid-correlation mid-weight robust to outliers
#' @importFrom Rsamtools indexFa index
#' @importFrom Rsamtools scanFaIndex FaFile scanFa path
#' @importFrom GenomeInfoDb seqnames
#' @importFrom DBI dbConnect dbDriver dbWriteTable dbGetQuery dbDisconnect
#' @importFrom S4Vectors DataFrame
#' @export
#' @return a SQLlite Data base 
repeatTablesFromWGCNA<-function(kexp,lnames,useBiCor=TRUE,verbose=TRUE,dbname="wgcnaDBLite",annotate=TRUE,version="1.0.0" ){
  
 ##FIX ME: create accessors for this class 
 ## FIX ME : adding this database will require a formal methods section
   tnxomes<-transcriptomes(kexp)
   ##supporting only ENSEMBL
   tnxomes<-strsplit(tnxomes,",")
   Ensmbl<-unlist(tnxomes)[grep("Ens",unlist(tnxomes))]
  
 if(verbose) cat("Extracting Networking data...\n")
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

  weight<-as.data.frame(datTraits)

  ##FIX ME: add organism
  if(all(grepl("^ENSG",colnames(datTraits)))==TRUE){
   organism<-"Homo.sapiens"
   packageName<-gsub("EnsDbLite","wgcnaDbLite",Ensmbl)
  } else if(all(grepl("^ENSMUSG",colnames(datTraits)))==TRUE){
  organism<-"Mus.musculus"
  packageName<-gsub("EnsDbLite","wgcanDbLite",Ensmbl)
  }





 if(useBiCor==FALSE){
  ##finds correlation per gene in modules
 if(verbose) cat("Calculating Module Correlations with EigenGenes\n")
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ##pvalue per gene per module
  if(verbose) cat("Calculating student-t pvalue ...\n")
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  if(verbose) cat("Calculating Cross Correlation of genes~traits...\n")
  geneTraitCor = as.data.frame(cor(datExpr, weight, use = "p"));
  if(verbose) cat("Calculating student-t pvalues of gene~trait correlations...\n")
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
  if(verbose) cat("Calculating Fisher-Exact p.values...\n")  
  geneTraitFisherPvalue = as.data.frame(corPvalueFisher(as.matrix(geneTraitCor), nSamples));
  colnames(geneTraitCor)<-paste("GCor.",colnames(weight),sep="");
  colnames(geneTraitPvalue)<-paste("p.GCor.",colnames(weight),sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf.GCr.",colnames(weight),sep="");
  } else {
  if(verbose) cat("Calculating Module bi-Correlations with EigenGenes...\n")
 geneModuleMembership = as.data.frame(bicor(datExpr, MEs, use = "all.obs"));
  if(verbose) cat("Calculating student-t p.values...\n")
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); ##pvalue per gene in each module
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  ##gene pvale per trait
  if(verbose) cat("Calculating bi-correlation of gene~trait...\n")
  geneTraitCor = as.data.frame(bicor(datExpr, weight, use = "all.obs"));
  if(verbose) cat("Calculating student-t p.values gene~trait...\n")
  geneTraitPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples));
  geneTraitFisherPvalue<-as.data.frame(corPvalueFisher(as.matrix(geneTraitCor),nSamples));
   colnames(geneTraitCor) = paste("GCor.", colnames(weight), sep="");
  colnames(geneTraitPvalue) = paste("p.GCor.", colnames(weight), sep="");
  colnames(geneTraitFisherPvalue)<-paste("pf.GCr.",colnames(weight),sep="");
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
  entrezID <- getEntrezIDs(annot, organism)
  entrez.match<-match(rownames(annot),names(entrezID))
  stopifnot(all(rownames(annot)==names(entrezID)[entrez.match]))
  annot$entrezid<-entrezID[entrez.match]
  symbolNames<-getSymbols(annot,organism)
  symbol.match<-match(rownames(annot),names(symbolNames))
  stopifnot(all(rownames(annot)==names(symbolNames)[symbol.match]))
  annot$hgnc_symbol<-symbolNames[symbol.match]
 # gxs$gene_name <- getSymbols(gxs, organism)[names(gxs)]
 }



 ##the goal is to have a SQL library instead of a large data frame geneTraitModuleDF
 ##the bwLabels is the key, and should be written to a table as well. this is the cornerstone to accessing moduleGenes and all associated data.

 if(verbose) cat("Creating Database...\n")
  
  dbname<-paste(dbname,"sqlite",sep=".")
  con<-dbConnect(dbDriver("SQLite"),dbname=dbname)
  if(verbose) cat("Extracting Module Color Names...")
   modulesNeeded<-unique(moduleColors)
   modulePvalue.Needed<-colnames(geneModuleDF)[grep("p.MM",colnames(geneModuleDF))]
  stopifnot(length(modulesNeeded)==length(modulePvalue.Needed))
  if(verbose) cat("Extracting Pheno-types of interest...")
  traitsNeeded<-colnames(datTraits)
  geneTraitCorNeeded<-paste0("GCor.",traitsNeeded)
  pTraitNeeded<-paste0("p.GCor.",traitsNeeded)
  pfTraitNeeded<-paste0("pf.GCr.",traitsNeeded) 
  stopifnot(all(geneTraitCorNeeded%in%colnames(geneModuleDF)))
  stopifnot(all(pTraitNeeded %in% colnames(geneModuleDF)))
  stopifnot(all(pfTraitNeeded%in%colnames(geneModuleDF)))
  
  geneTraitCor.id<-match(geneTraitCorNeeded,colnames(geneModuleDF))
  pTrait.id<-match(pTraitNeeded,colnames(geneModuleDF))
  pfTrait.id<-match(pfTraitNeeded,colnames(geneModuleDF))
  #color,GCor_TxBio,pvalue , pf.value,annotations
  ## FIX ME: add annotations
  annot.id<-match(rownames(geneModuleDF),annot$ensembl_gene_id)
  
  for(i in 1:length(tablesNeeded)){
  mod.id<-which(modulesNeeded[i]==substring(colnames(geneModuleDF),3 ))
  mod.p.id<-which(modulePvalue.Needed[i]==substring(colnames(geneModuleDF),5))
  colorTable<-geneModuleDF[,c(mod.id,mod.p.id,geneTraitCor.id,pTrait.id,pfTrait.id)]
  colorKey<-as.factor(moduleColors)
  colorTable<-cbind(colorTable,colorKey)
  colorTable<-cbind(colorTable,annot[annot.id,]) 
  stopifnot(all(rownames(colorTable)==rownames(geneModuleDF)))
  stopifnot(all(rownames(colorTable)==colorTable$ensembl_gene_id))
  
  ##write dbLite table  
  if(verbose) cat(paste0("Creating the database for ",tablesNeeded[i],"\n"))
  dbWriteTable(con,name=tablesNeeded[i],colorTable,overwrite=T,row.names=T)
  dbGetQuery(con,paste0("create index colorKey_idx","_",i," on ",tablesNeeded[i]," (colorKey);"))
   }
  ##FIX ME :  write out some meta data color module summaries, species .. 
  Metadata <- wgcnaDbLiteMetadata(kexp,
                               packageName=packageName,
                               organism=organism,
                               Ensmbl=Ensmbl)
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
#' @param organism  what kind of organism these genes are from
#' 
#' @return  entrez_id values for the genes, where found
#'
#'
getEntrezIDs <- function(gxs, organism) { # {{{
  org <- getOrgDetails(organism)
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
#' add symbols for Ensembl genes
#' @return  symbols for the genes, where found 
#'
#' 
getSymbols <- function(gxs, organism) { # {{{
  org <- getOrgDetails(organism)
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
wgcnaDbLiteMetadata <- function(kexp,packageName,organism,Ensmbl) { # {{{

  
  Ensmbl<-gsub(" ","",Ensmbl)
  org <- getOrgDetails(organism)
  sourceFile <- Ensmbl
  kversion<-kallistoVersion(kexp)
  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "wgcnaDbLite")
  MetaData[3,] <- c("type_of_gene_id", "Ensembl Gene ID")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", organism)
  MetaData[7,] <- c("kallistoVersion", kversion)
  MetaData[8,] <- c("sourceFile", sourceFile)
  return(MetaData)

} # }}}

