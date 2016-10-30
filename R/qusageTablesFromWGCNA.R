#' @title qusage enrichme modules returned from WGCNA.  the goal is to then select a MSigDB gmt and take the module.tx_biotype of interest and create a corresponding qusage enrichment database.
#' @param kexp the kexp should have pHSC, LSC and blast.  
#' @param verbose boolean for std out
#' @param dbname this is the db name for wgcnaDbLite with correlation db
#' @param version character, for the metadata information
#' @param Module.color this is a color needed to grab the trait names 
#' @param MsigDB which db to use
#' @param how either cpm or tpm
#' @param species either Homo.sapiesn or Mus.musculus
#' @param comparison1 pHSC, or charactoer for 2 group copmarison. these arguments are input into qusageRun and must be split using "_" with the factor of interest as the first leading term.
#' @param controls LSC, or character for control group. these parameters are input into qusageRun with the factor of interest as the leading term.
#' @param comparison2 Blast, or character for 2nd pairwise. this is required, if you only have a two group then let comparison1=comparison2.  must be named using "_" with the factor of interest as the leading term.
#' @param paired for patient paired data
#' @param batchNormalize boolean, if true then the samples are from different batches and the edgeR batchCorrectEffect will be called for cpm
#' @param batchVector this is the batch vector. the length of the batch vector must match the number of columns. there are not any checks to validate the batch effect factors, user dependent input
#' @param comparisonNumber integer 1 or 2.  this integer reflects how many pairwise comparisons needed.  for time series analysis, it is sometimes desired to test for enrichment at N different time points,  intial-comparison1 , and intial-comparison2. qusage will be called to compare the comparisonNumber pairs 
#' @import edgeR
#' @import TxDbLite
#' @export
#' @return a qusageDbLite db
qusageTablesFromWGCNA<-function(kexp,verbose=TRUE,dbname="wgcnaDBLite.sqlite",version="1.0.0",Module.color="brown",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),how=c("cpm","tpm"),species=c("Homo.sapiens","Mus.musculus"),comparison1="pHSC",comparison2="Blast",controls="LSC",paired=TRUE,batchNormalize=FALSE,batchVector=NULL,comparisonNumber=1 ){

##task: this will take a full kexp and first normalize then collapse by gene_name then log2 transform or tpm normalize and pass into qusageRun.R to handle the pre-proccessing for qusage call.  the output should be module.biotype.enrich specific data.  then write to a db. 
 ##important note: it is tempting to merely use collapseBundles(kexp,"gene_name") and pipe into qusage HOWEVER QUSAGE REQUIRES NORMALIZED LOG2 XR. collapseBundles will only use the raw counts(kexp) bundle gene_name counts.  here we *****MUST******* use tmm normalized or Tpm calls.
##ALSO NOTE: you can first collapse by gene_id normalize w.r.t gene_id and then collapse by gene name, but that is not that same as collapsing by gene_id, querying a module, then collapsing by gene name and then normalizing by gene_name.  we see better results for the latter


   tnxomes<-transcriptomes(kexp)
   ##supporting only ENSEMBL
   tnxomes<-strsplit(tnxomes,",")
   Ensmbl<-unlist(tnxomes)[grep("Ens",unlist(tnxomes))]
   how<-match.arg(how,c("cpm","tpm"))
   geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))

##### the modules are in terms of ENSG IDs and we are going to build a qusage database of the modules, so we must collapse counts by ENSG IDs, subset the module, collapse module genes by gene_name, and then normalize, log2XR and then annotate to get hgnc to pipe into qusage.  
#####
  
 
    if(how=="cpm"){
  ##TMM normalize
  full_counts<-collapseBundles(kexp,"gene_id")
 # dge<-DGEList(counts=counts)
 # dge<-calcNormFactors(dge)
 # expr<-cpm(dge,log=FALSE)
 # counts<-log2(1+expr) ##enirhcment on log2 is required
  } else{
  ##TPM Normalize
  full_counts<-collapseTpm(kexp,"gene_id")
   full_counts<-log2(1+counts)
  }

  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
   if(species=="Homo.sapiens"){
   packageName<-gsub("EnsDbLite","wgcnaDbLite",Ensmbl)
  } else if(species=="Mus.musculus" ){
  packageName<-gsub("EnsDbLite","wgcanDbLite",Ensmbl)
  }

 
 ## take the database and query the geneIDs of interest, and use the kexp on those geneIDs
   
 allTables<-dbListTables(dbconn(wgcnaDbLite(dbname)))
  ##the tables have metadata, must avoid meta table
 allcolors<-allTables[!grepl("metadata",allTables)] 
 allcolors<-allcolors[!grepl("go",allcolors)] ##filter out goEnrich table
 allTraits<-colnames(modulesBy(wgcnaDbLite(dbname),Module.color=Module.color)) 
 allTraits<-allTraits[grepl("^GCor_",allTraits)]
 allTraits<-sapply(strsplit(allTraits,"_"),function(x) x[2])
  if(verbose) cat("Creating Enrichment database...\n")
  qusage.dbname<-paste0("qusageDbLite.",how,".sqlite")
  quscon<-dbConnect(dbDriver("SQLite"),dbname=qusage.dbname)
  
## take the module genes unfiltered --> qusage (modulesBy)
 for(i in 1:length(allcolors)){
 ######START FOR LOOP HERE allcolors_i 
   wgcna.color<-modulesBy(wgcnaDbLite(dbname),p.value=2,Module.color=allcolors[i])
   color.id<-match(wgcna.color$row_names,rownames(full_counts))
   stopifnot(all(rownames(full_counts[color.id,] )==wgcna.color$row_names)) ##check
   ##subset module by ENSGID
  module.counts<-full_counts[color.id,]
  module.counts<-module.counts[,which(colSums(module.counts)>1)]
  ## first collapse by gene name and then normalize.  as opposed to normalizing gene_ids and then collapsing the modules.
   ##module kexp
  module.kexp<-kexp[rowRanges(kexp)$gene_id%in%rownames(module.counts),]
   if(how=="cpm"){
   counts<-collapseBundles(module.kexp,"gene_name")
   counts<-counts[,which(colSums(counts)>1)]
   dge<-DGEList(counts=counts)
   dge<-calcNormFactors(dge)
   expr<-cpm(dge,normalized.lib.sizes=TRUE,log=FALSE)
   if(batchNormalize==TRUE){
     stopifnot(is.null(batchVector)==FALSE)
     stopifnot(length(batchVector)==ncol(kexp))##the batchVector nomenclature must match the length of columns.
     ##takes TMM normalized and batch corrects
     batch.cpm<-removeBatchEffect(expr,batch=batchVector)
      }
   if(batchNormalize==FALSE){ 
  counts<-log2(1+expr) ##enirhcment on log2 is required
  }else if(batchNormalize==TRUE){
   counts<-log2(1+batch.cpm)
    ### filter NaNs
    id<-which(is.na(counts))
    counts[id]<-1.0
    }
  } else {
  counts<-collapseTpm(module.kexp,"gene_name")
  counts<-log2(1+counts)
   }  
   ##split the count data by stage
   module.counts<-counts
   comparison1.id<-grep(paste0(comparison1,"_"),colnames(module.counts))
   controls.id<-grep(paste0(controls,"_"),colnames(module.counts))
   comparison2.id<-grep(paste0(comparison2,"_"),colnames(module.counts))
   cnts_phLS<-module.counts[,c(comparison1.id,controls.id)]
   cnts_blsLS<-module.counts[,c(comparison2.id,controls.id)]

    ##call qusage for each stage
    qusage_run1<-qusageRun(cnts_mt=cnts_phLS,MsigDB=MsigDB,comparison=comparison1,control=controls,module=allcolors[i],paired=paired)
     if(comparisonNumber==2){
    qusage_run2<-qusageRun(cnts_mt=cnts_blsLS,MsigDB=MsigDB,comparison=comparison2,control=controls,module=allcolors[i],paired=paired)
    }else if(comparisonNumber==1){
     qusage_run2<-qusage_run1
     #if comparisonNumber1 is true then we check and use only comparison1-controls and copy that into a dummy run2.
    }else{
     cat("currently support only 2 paired time points of stage.\n")
    }
      if(nrow(qusage_run1[!is.na(qusage_run1$p.Value),])>1){
      qusage_run1<-data.frame(qusage_run1[!is.na(qusage_run1$p.Value),],
                             colorKey=allcolors[i],
                             bioKey=allcolors[i],
                             contrastKey=tolower(comparison1))
       }else{
         qusage_run1<-data.frame(pathway.name="NA",
                                    log.fold.change="NA",
                                    p.Value="NA",
                                    FDR="NA",
                                    colorKey=allcolors[i],
                                    bioKey=allcolors[i],
                                    contrastKey=tolower(comparison1))
       } 
      if(nrow(qusage_run2[!is.na(qusage_run2$p.Value),])>1){                             qusage_run2<-data.frame(qusage_run2[!is.na(qusage_run2$p.Value),] ,
                              colorKey=allcolors[i], 
                              bioKey=allcolors[i],
                              contrastKey=tolower(comparison2))

       }else{
         qusage_run2<-data.frame(pathway.name="NA",
                                    log.fold.change="NA",
                                    p.Value="NA",
                                    FDR="NA",
                                    colorKey=allcolors[i],
                                    bioKey=allcolors[i],
                                    contrastKey=tolower(comparison2))

       }
   ##entire module universe of enrichment.
          qusage_module<-rbind(qusage_run1,qusage_run2) 
    for( j in 1:length(allTraits)){
   ########SECOND FOR LOOP ACROSS ALL TRAITS HERE #########
   ###here we filter the most significant correlated module ~ trait (driver) associations and run enrichment on these
  if(verbose) cat("Extracting Trait Driver Enrichments...\n")  
  ##grabs the wgcnaDbLite db 
 wgcna.color.trait<-traitsBy(wgcnaDbLite(dbname),p.value=0.05,Module.color=allcolors[i],trait=allTraits[j])
   
  if(nrow(wgcna.color.trait[!is.na(wgcna.color.trait$hgnc_symbol),]  )<=8){
    qusage_module.trait<-data.frame(pathway.name="NA",
                                    log.fold.change="NA",
                                    p.Value="NA",
                                    FDR="NA",
                                    colorKey=allcolors[i],
                                    bioKey=allTraits[j],
                                    contrastKey="NA")
    qusage_module<-rbind(qusage_module,qusage_module.trait)
   next 
   } 
   trait.color.id<-match(wgcna.color.trait$row_names,rownames(full_counts))
   stopifnot(rownames(full_counts[trait.color.id,])==wgcna.color.trait$row_names)
   module.trait.counts<-full_counts[trait.color.id,]
    ## collapse by gene name and then normalize moduleKexp.  as opposed to normalizing first and then mapping to gene names.
  module.trait.kexp<-kexp[rowRanges(kexp)$gene_id%in%rownames(module.trait.counts),]
   if(how=="cpm"){
   module.trait.counts<-collapseBundles(module.trait.kexp,"gene_name")
   module.trait.counts<-module.trait.counts[,which(colSums(module.trait.counts)>1)]
   if(nrow(module.trait.counts)<=8){
   next
   }
   trait.dge<-DGEList(counts=module.trait.counts)
   trait.dge<-calcNormFactors(trait.dge)
   trait.expr<-cpm(trait.dge,log=FALSE)
   module.trait.counts<-log2(1+trait.expr) ##enirhcment on log2 is required
    } else {
  module.trait.counts<-collapseTpm(module.trait.kexp,"gene_name")
   if(nrow(module.trait.counts)<=8){
    next
   }
  module.trait.counts<-log2(1+module.trait.counts)
   }
   ##split the count data by stage
   comparison1.trait.id<-grep(paste0(comparison1,"_"),colnames(module.trait.counts))
   controls.trait.id<-grep(paste0(controls,"_"),colnames(module.trait.counts))
   comparison2.trait.id<-grep(paste0(comparison2,"_"),colnames(module.trait.counts))
   trait.cnts_phLS<-module.trait.counts[,c(comparison1.trait.id,controls.trait.id)]
   trait.cnts_blsLS<-module.trait.counts[,c(comparison2.trait.id,controls.trait.id)]
  stopifnot(nrow(trait.cnts_phLS)==nrow(module.trait.counts))
  stopifnot(nrow(trait.cnts_blsLS)==nrow(module.trait.counts))
    ##call qusage for each stage
    qusage_trait.run1<-qusageRun(cnts_mt=trait.cnts_phLS,MsigDB=MsigDB,comparison=comparison1,control=controls,module=allcolors[i],paired=paired)
    qusage_trait.run2<-qusageRun(cnts_mt=trait.cnts_blsLS,MsigDB=MsigDB,comparison=comparison2,control=controls,module=allcolors[i],paired=paired)
     
     if(nrow(qusage_trait.run1[!is.na(qusage_trait.run1$p.Value),])>1){
     qusage_trait.run1<-data.frame(qusage_trait.run1[!is.na(qusage_trait.run1$p.Value),],
                             colorKey=allcolors[i],
                             bioKey=allTraits[j],
                             contrastKey=tolower(comparison1))
     }else{
        qusage_trait.run1<-data.frame(pathway.name="NA",
                                    log.fold.change="NA",
                                    p.Value="NA",
                                    FDR="NA",
                                    colorKey=allcolors[i],
                                    bioKey=allTraits[j],
                                    contrastKey=tolower(comparison1))

     }
     if(nrow(qusage_trait.run2[!is.na(qusage_trait.run2$p.Value),])>1){
       qusage_trait.run2<-data.frame(qusage_trait.run2[!is.na(qusage_trait.run2$p.Value),],
                              colorKey=allcolors[i],
                              bioKey=allTraits[j],
                              contrastKey=tolower(comparison2))

     }else{
     qusage_trait.run2<-data.frame(pathway.name="NA",
                                    log.fold.change="NA",
                                    p.Value="NA",
                                    FDR="NA",
                                    colorKey=allcolors[i],
                                    bioKey=allTraits[j],
                                    contrastKey=tolower(comparison2))

       }
   
      qusage_module.trait<-rbind(qusage_trait.run1,qusage_trait.run2)
  
      qusage_module<-rbind(qusage_module,qusage_module.trait)

  } ##end of j allTrait loop
  if(verbose) cat("Writing Module Enrichment Database...\n")

  colnames(qusage_module)<-c("pathway_name","logFC","pvalue","FDR","colorKey","bioKey","contrastKey")
   dbWriteTable(quscon,name=allcolors[i],qusage_module,overwrite=T,row.names=T)
  
 }## END ALL i allcolors loop
  if(verbose) cat("Writing enrichment of entire kexp...\n")
  
 ##FIX ME: call for the ***whole*** pHSC-LSC and Blast-LSC as separate tables
 ## take the database and query the geneIDs of interest, and use the kexp on those geneIDs
  if(how=="cpm"){
  ##TMM normalize
  full_counts<-collapseBundles(kexp,"gene_name")
  dge<-DGEList(counts=full_counts)
  dge<-calcNormFactors(dge)
  expr<-cpm(dge,log=FALSE)
  full_counts<-log2(1+expr) ##enirhcment on log2 is required
  } else{
  ##TPM Normalize
  full_counts<-collapseTpm(kexp,"gene_id")
   full_counts<-log2(1+full_counts)
  }
 ##split by stage
   full_module.counts<-full_counts
   pHSC.id<-grep(paste0(comparison1,"_"),colnames(full_module.counts))
   LSC.id<-grep(paste0(controls,"_"),colnames(full_module.counts))
   Blast.id<-grep(paste0(comparison2,"_"),colnames(full_module.counts))
   full_cnts_phLS<-full_module.counts[,c(pHSC.id,LSC.id)]
   full_cnts_blsLS<-full_module.counts[,c(Blast.id,LSC.id)]

    ##call qusage for each stage
    qusage_full1<-qusageRun(cnts_mt=full_cnts_phLS,MsigDB=MsigDB,comparison=comparison1,control=controls,module="fullKexp",paired=paired)
    qusage_full2<-qusageRun(cnts_mt=full_cnts_blsLS,MsigDB=MsigDB,comparison=comparison2,control=controls,module="fullKexp",paired=paired)
      qusage_full1<-data.frame(qusage_full1,
                             colorKey="kexp",
                             bioKey="kexp",
                             contrastKey=tolower(comparison1))
    
        qusage_full2<-data.frame(qusage_full2,
                              colorKey="kexp",
                              bioKey="kexp",
                              contrastKey=tolower(comparison2))
  full_qusage<-rbind(qusage_full1,qusage_full2)
  dbWriteTable(quscon,name="kexp",full_qusage,overwrite=T,row.names=T)


 ##writeTable
 

 if(verbose) cat("Writing Metadata table...\n")
   Metadata <- qusageDbLiteMetadata(kexp,
                               packageName=packageName,
                               species=species,
                               Ensmbl=Ensmbl,
                               how=how)
 dbWriteTable(quscon, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)


  cat("done.\n")
  dbDisconnect(quscon)
  return(qusage.dbname)
} ##main




#' @describeIn qusageTablesFromWGCNA 
#' 
#' create metadata for an qusageDbLite instance
#'
#' @param packageName   the name of the annotation package to be built 
#' @param genomeVersion name of genome assembly for coordinates, e.g. "GRCh38"
#' @param sourceFile    name of FASTA file(s) whence it was built, as a string
#' 
#' @return a data.frame of metadata suitable for cramming into the database
#' 
#' @export
qusageDbLiteMetadata <- function(kexp,packageName,species=NULL,Ensmbl,how=how) { # {{{


  Ensmbl<-gsub(" ","",Ensmbl)
  org <- getOrgDetails(species)
  sourceFile <- Ensmbl
  kversion<-kallistoVersion(kexp)
  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "qusageDbLite")
  MetaData[3,] <- c("type_of_gene_id", "Ensembl Gene ID")
  MetaData[4,] <- c("created_by", how)
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", species)
  MetaData[7,] <- c("qusage_Version", paste("qusage",packageVersion("qusage")))
  MetaData[8,] <- c("sourceFile", sourceFile)
  return(MetaData)

} # }}}


