#' @title qusage enrichme modules returned from WGCNA.  the goal is to then select a MSigDB gmt and take the module.tx_biotype of interest and create a corresponding qusage enrichment database.
#' @param kexp the kexp should have pHSC, LSC and blast.  
#' @param verbose boolean for std out
#' @param dbname this is the db name for wgcnaDbLite with correlation db
#' @param version character, for the metadata information
#' @param Module.color this is a color needed to grab the trait names 
#' @param MsigDB which db to use
#' @param how either cpm or tpm
#' @param species either Homo.sapiesn or Mus.musculus
#' @import edgeR
#' @import TxDbLite
#' @export
#' @return a qusageDbLite db
qusageTablesFromWGCNA<-function(kexp,verbose=TRUE,dbname="wgcnaDBLite.sqlite",version="1.0.0",Module.color="brown",MsigDB=c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"),how=c("cpm","tpm"),species=c("Homo.sapiens","Mus.musculus") ){

##task: this will take a full kexp and collapse by gene_name then  tmm normalize and log2 transform or tpm normalize and pass into qusageRun.R to handle the pre-proccessing for qusage call.  the output should be module.biotype.enrich specific data.  then write to a db. 
 ##important note: it is tempting to merely use collapseBundles(kexp,"gene_name") and pipe into qusage HOWEVER QUSAGE REQUIRES NORMALIZED LOG2 XR. collapseBundles will only use the raw counts(kexp) bundle gene_name counts.  here we *****MUST******* use tmm normalized or Tpm calls.
##ALSO NOTE: you can first normalize and then collapse by gene name, but that is not that same as collapsing by gene name and then normalizing.  we see better results for the latter.  


   tnxomes<-transcriptomes(kexp)
   ##supporting only ENSEMBL
   tnxomes<-strsplit(tnxomes,",")
   Ensmbl<-unlist(tnxomes)[grep("Ens",unlist(tnxomes))]
   how<-match.arg(how,c("cpm","tpm"))
   geneSet<-match.arg(MsigDB,c("c1.all.v5.1.symbols.gmt","c2.all.v5.1.symbols.gmt","c4.all.v5.1.symbols.gmt","c5.all.v5.1.symbols.gmt","c6.all.v5.1.symbols.gmt","c7.all.v5.1.symbols.gmt","h.all.v5.1.symbols.gmt"))

##### the modules are in terms of ENSG IDs and we are going to build a qusage database of the modules, so we must collapse counts by ENSG IDs,normalize, log2XR and then annotate to get hgnc to pipe into qusage.  
#####
  
  ##Lazy Species detection
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

  ##FIX ME: call for the whole pHSC-LSC and Blast-LSC as separate tables
 ## take the database and query the geneIDs of interest, and use the kexp on those geneIDs
   
 allTables<-dbListTables(dbconn(wgcnaDbLite(dbname)))
  ##the tables have metadata, must avoid this table
 allcolors<-allTables[!grepl("metadata",allTables)]
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
module.counts<-full_counts[color.id,]
  ##we found better results for first collapse by gene name and then normalizing.  as opposed to normalizing first and then mapping to gene names.
  module.kexp<-kexp[rowRanges(kexp)$gene_id%in%rownames(module.counts),]
   if(how=="cpm"){
   counts<-collapseBundles(module.kexp,"gene_name")
   dge<-DGEList(counts=counts)
   dge<-calcNormFactors(dge)
   expr<-cpm(dge,log=FALSE)
   counts<-log2(1+expr) ##enirhcment on log2 is required
    } else {
  counts<-collapseTpm(module.kexp,"gene_name")
  counts<-log2(1+counts)
   }  
   ##split the count data by stage
   module.counts<-counts
   pHSC.id<-grep("pHSC_",colnames(module.counts))
   LSC.id<-grep("LSC_",colnames(module.counts))
   Blast.id<-grep("Blast_",colnames(module.counts))
   cnts_phLS<-module.counts[,c(pHSC.id,LSC.id)]
   cnts_blsLS<-module.counts[,c(Blast.id,LSC.id)]

    ##call qusage for each stage
    qusage_run1<-qusageRun(cnts_mt=cnts_phLS,MsigDB=MsigDB,comparison="pHSC",control="LSC",module=allcolors[i])
    qusage_run2<-qusageRun(cnts_mt=cnts_blsLS,MsigDB=MsigDB,comparison="Blast",control="LSC",module=allcolors[i])
      qusage_run1<-data.frame(qusage_run1[!is.na(qusage_run1$p.Value),],
                             colorKey=allcolors[i],
                             bioKey=allcolors[i],
                             contrastKey="phsc")
                                
       qusage_run2<-data.frame(qusage_run2[!is.na(qusage_run2$p.Value),] ,
                              colorKey=allcolors[i], 
                              bioKey=allcolors[i],
                              contrastKey="blast")
   ##entire module universe of enrichment.
          qusage_module<-rbind(qusage_run1,qusage_run2) 
    for( j in 1:length(allTraits)){
   ########SECOND FOR LOOP ACROSS ALL TRAITS HERE #########
   ###here we filter the most significant correlated module ~ trait driver and run enrichment on these
  if(verbose) cat("Extracting Trait Driver Enrichments...\n")  
  ##grabs the wgcnaDbLite db 
 wgcna.color.trait<-traitsBy(wgcnaDbLite(dbname),p.value=0.05,Module.color=allcolors[i],trait=allTraits[j]) 
  if(nrow(wgcna.color.trait)<=8){
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
    ##we found better results for first collapse by gene name and then normalizing.  as opposed to normalizing first and then mapping to gene names.
  module.trait.kexp<-kexp[rowRanges(kexp)$gene_id%in%rownames(module.trait.counts),]
   if(how=="cpm"){
   module.trait.counts<-collapseBundles(module.trait.kexp,"gene_name")
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
   pHSC.trait.id<-grep("pHSC_",colnames(module.trait.counts))
   LSC.trait.id<-grep("LSC_",colnames(module.trait.counts))
   Blast.trait.id<-grep("Blast_",colnames(module.trait.counts))
   trait.cnts_phLS<-module.trait.counts[,c(pHSC.trait.id,LSC.trait.id)]
   trait.cnts_blsLS<-module.trait.counts[,c(Blast.trait.id,LSC.trait.id)]
  stopifnot(nrow(trait.cnts_phLS)==nrow(module.trait.counts))
  stopifnot(nrow(trait.cnts_blsLS)==nrow(module.trait.counts))
    ##call qusage for each stage
    qusage_trait.run1<-qusageRun(cnts_mt=trait.cnts_phLS,MsigDB=MsigDB,comparison="pHSC",control="LSC",module=allcolors[i])
    qusage_trait.run2<-qusageRun(cnts_mt=trait.cnts_blsLS,MsigDB=MsigDB,comparison="Blast",control="LSC",module=allcolors[i])
      qusage_trait.run1<-data.frame(qusage_trait.run1,
                             colorKey=allcolors[i],
                             bioKey=allTraits[j],
                             contrastKey="phsc")

       qusage_trait.run2<-data.frame(qusage_trait.run2,
                              colorKey=allcolors[i],
                              bioKey=allTraits[j],
                              contrastKey="blast")
      qusage_module.trait<-rbind(qusage_trait.run1,qusage_trait.run2)
  
      qusage_module<-rbind(qusage_module,qusage_module.trait)

  } ##end of j allTrait loop
  if(verbose) cat("Writing Module Enrichment Database...\n")

  colnames(qusage_module)<-c("pathway_name","logFC","pvalue","FDR","colorKey","bioKey","contrastKey")
   dbWriteTable(quscon,name=allcolors[i],qusage_module,overwrite=T,row.names=T)
  
 }## END ALL i allcolors loop
  if(verbose) cat("Writing enrichment of entire kexp...\n")
  
 ##FIX ME: call for the whole pHSC-LSC and Blast-LSC as separate tables
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
   pHSC.id<-grep("pHSC_",colnames(full_module.counts))
   LSC.id<-grep("LSC_",colnames(full_module.counts))
   Blast.id<-grep("Blast_",colnames(full_module.counts))
   full_cnts_phLS<-full_module.counts[,c(pHSC.id,LSC.id)]
   full_cnts_blsLS<-full_module.counts[,c(Blast.id,LSC.id)]

    ##call qusage for each stage
    qusage_full1<-qusageRun(cnts_mt=full_cnts_phLS,MsigDB=MsigDB,comparison="pHSC",control="LSC",module="fullKexp")
    qusage_full2<-qusageRun(cnts_mt=full_cnts_blsLS,MsigDB=MsigDB,comparison="Blast",control="LSC",module="fullKexp")
      qusage_full1<-data.frame(qusage_full1,
                             colorKey="kexp",
                             bioKey="kexp",
                             contrastKey="phsc")
    
        qusage_full2<-data.frame(qusage_full2,
                              colorKey="kexp",
                              bioKey="kexp",
                              contrastKey="blast")
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


