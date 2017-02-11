#' @rdname wgcnaDbLite-class
#' @export
setMethod("dbconn", "wgcnaDbLite", function(x) return(x@con))

#' @rdname wgcnaDbLite-class
#' @param object wgcnaDbLite SQL-lite annotation database
#' @import ensembldb
#' @export
setMethod("show", "wgcnaDbLite", function(object) { # {{{
  if(is.null(object@con)) stop(paste("Invalid", class(object), "instance!"))
  info <- metadata(object)
  cat(class(object), ":\n")
  catrow <- function(x) cat(paste0("|", x["name"], ": ", x["value"], "\n"))
  for(i in 1:nrow(info)) catrow(info[i,])
}) # }}}

#' @rdname wgcnaDbLite-class
#' @param x wgcnaDblite SQL-lite database instance
#' @param ... additional parameters releated to annotation database
#' @export
setMethod("metadata", "wgcnaDbLite", function(x, ...) { # {{{
  md <- dbGetQuery(dbconn(x), "select * from metadata")
  rownames(md) <- md$name
  return(md)
}) # }}}


setGeneric("modulesBy",function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL,...) standardGeneric("modulesBy")) 
 
  
#' @rdname wgcnaDbLite-class
#' @description this will return all the gene cross correlations for a given module and all traits in a given row if trait param is NULL. or if traits are specified will return the cross correlations (significant) between a desired module and trait.
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param trait  this is the trait of interest 
#' @param p.value p.value significance threshold
#' @export
setMethod("modulesBy", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL) {
  p.value.type<-match.arg(p.value.type,c("studentT","fisher"))
  if(p.value.type=="studentT"){
  p.prefix<-"p_GCor_"
  alt.p.prefix<-"pf_GCr_"
  } else{
  p.prefix<-"pf_GCr_"
  alt.p.prefix<-"p_GCor_"
  }
  ##FIX ME: safety check for the input color using dbListTables,  signature
  allcolors<-dbListTables(dbconn(x))
  stopifnot(Module.color%in%allcolors ==TRUE)
  if(is.null(trait)==TRUE){
  #trait was not specified so return all traits 
  sql<-paste0("select * from ",Module.color," where colorKey='",Module.color,"'")
   res<-as.data.frame(dbGetQuery(dbconn(x),sql),
                    stringsAsFactors=FALSE)
   ##only return one specified p.value type 
   id<-!grepl(alt.p.prefix,colnames(res))
   res<-res[,id]
  } else{
   #traits was specified so return specific
  gcor<-paste0("GCor_",trait)
  pval.Trait<-paste0(p.prefix,trait)
  sql<-paste0("select row_names, ","MM",Module.color,", ",gcor,", ",pval.Trait,", colorKey, gene_id, entrezid, hgnc_symbol from ",Module.color," where colorKey='",Module.color,"'")
  res<-as.data.frame(dbGetQuery(dbconn(x),sql))
     if(is.null(p.value)==FALSE){
    p.col<-paste0(p.prefix,trait)
    res<-res[which(res[,grepl(p.col,colnames(res))]<=p.value),]
   } 
   }
    return(res)
  
 })



setGeneric("traitsBy",function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL,...) standardGeneric("traitsBy"))



#' @rdname wgcnaDbLite-class
#' @description traitsBy will return significant Gene cross correlations  between a trait and module of interest. this is synonomous to modulesBy with a specifid p.value and trait.  
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param trait  this is the trait of interest 
#' @param p.value p.value significance threshold
#' @export
setMethod("traitsBy", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL) {
  p.value.type<-match.arg(p.value.type,c("studentT","fisher"))
  if(p.value.type=="studentT"){
  p.prefix<-"p_GCor_"
  alt.p.prefix<-"pf_GCr_"
  } else{
  p.prefix<-"pf_GCr_"
  alt.p.prefix<-"p_GCor_"
  }
  ##FIX ME: safety check for the input color using dbListTables,  signature
  allcolors<-dbListTables(dbconn(x))
  stopifnot(Module.color%in%allcolors ==TRUE)
  
  #drivers are defined as the most significantly cross correlated genes
  gcor<-paste0("`","GCor_",trait,"`")
  pval.Trait<-paste0("`",p.prefix,trait,"`")
  sql<-paste0("select row_names, ","MM",Module.color,", ",gcor,", ",pval.Trait,", colorKey, gene_id, entrezid, hgnc_symbol from ",Module.color," where colorKey='",Module.color,"'")
  res<-as.data.frame(dbGetQuery(dbconn(x),sql))
     if(is.null(p.value)==TRUE){
    p.col<-paste0(p.prefix,trait)
    res<-res[which(res[,grepl(p.col,colnames(res))]<=0.05),]
   } else{
     p.col<-paste0(p.prefix,trait)
    res<-res[which(res[,grepl(p.col,colnames(res))]<=p.value),]
    }
  
    return(res)

 })

setGeneric("drivers",function(x,Module.color=NULL,p.value.type="studentT",p.value=NULL,traitConsensus=NULL,species=NULL,...) standardGeneric("drivers"))


#' @rdname wgcnaDbLite-class
#' @description we define a 'driver' of a module~trait association as the set of all genes shared amongst correlated sets with a significance. also for the intersected genes across traits, the correlation scores for each gene are averages, along with the p.values so that the intersection 'driver' list will have an average correlation score, and an average p.value for a driver association across sets of traits.
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param p.value p.value significance threshold
#' @export
setMethod("drivers", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",p.value=0.05,traitConsensus=c("Alu","ERVK","ERV3","ERVL","LTR Retrotransposon"),species="Homo.sapiens") {

  p.value.type<-match.arg(p.value.type,c("studentT","fisher"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))

  if(p.value.type=="studentT"){
  p.prefix<-"p_GCor_"
  alt.p.prefix<-"pf_GCr_"
  } else{
  p.prefix<-"pf_GCr_"
  alt.p.prefix<-"p_GCor_"
  }
  ##FIX ME: safety check for the input color using dbListTables,  signature
  allcolors<-dbListTables(dbconn(x))
  stopifnot(Module.color%in%allcolors ==TRUE)
  
  #trait was not specified so return all traits 
  sql<-paste0("select * from ",Module.color," where colorKey='",Module.color,"'")

   res<-as.data.frame(dbGetQuery(dbconn(x),sql),
                    stringsAsFactors=FALSE)
   ##only return one specified p.value type 
   id<-!grepl(paste0("^",alt.p.prefix,"*"),colnames(res))
   #select one module and every trait (column)
   res<-res[,id]
 
    
   ###  find the IDs of columns
   traits.id<-match(traitConsensus,substring(colnames(res),6))
   stopifnot(length(traits.id)>0)
   module.id<-grep(paste0("MM",Module.color),colnames(res))
   rownames.id<-grep("row_names",colnames(res))
   #select the first 
   seed<-1
   for(i in 1:length(traits.id)){
   col.id<-which(substring(colnames(res)[traits.id[i]],6)==substring(colnames(res),8))
    if(i==seed){
    #we grab the first trait of interest, and filter by p.value, and must ensure that it is non-empty
    out<-res[,c(rownames.id,module.id,traits.id[i],col.id)] 
    p.col<-paste0(p.prefix,substring(colnames(res)[traits.id[i]],6))
    out<-out[which(out[,grepl(p.col,colnames(out))]<=p.value),]
    out.GCor.id<-grep("^GCor_",colnames(out))
    out.p.id<-grep("^p_GCor_",colnames(out))
 
   if(nrow(out)<=0){
     cat(paste0("did not detect significant correlations for ",colnames(res)[traits.id[i]],"\n"))
      seed<-seed+1
      next 
     }# first i
   } else if(i>seed){
    ##take the next iteration intersect, and average and repeat across all traits
    out2<-res[,c(rownames.id,module.id,traits.id[i],col.id)] 
    p.col2<-paste0(p.prefix,substring(colnames(res)[traits.id[i]],6))
    out2<-out2[which(out2[,grepl(p.col2,colnames(out2))]<=p.value),]
    out2.GCor.id<-grep("^GCor_",colnames(out2))
    out2.p.id<-grep("^p_GCor_",colnames(out2))
    interTraits<-intersect(out$row_names,out2$row_names)
    if(length(interTraits)<=0){
    next
    } else{
    gxs<-data.frame(gene_id=interTraits,stringsAsFactors=FALSE)
    rownames(gxs)<-gxs$gene_id
    out.id<-match(interTraits,out$row_names)
    out2.id<-match(interTraits,out2$row_names)
    stopifnot(out[out.id,]$row_names==out2[out2.id,]$row_names)
     ##average GCor and average p.value
    drivers<-data.frame(row_names=interTraits,
               avg_GCor=(out[out.id,out.GCor.id]+out2[out2.id,out2.GCor.id])/2,
                avg_p.value=(out[out.id,out.p.id]+out2[out2.id,out2.p.id])/2)
    ###important update objects!!!!!
    out<-drivers
    out.GCor.id<-grep("^avg_GCor",colnames(out))
    out.p.id<-grep("avg_p.value",colnames(out))

   }
  }
  }#for
  
 
    return(drivers)

 })


setGeneric("pathways",function(x,Module.color=NULL,tx.Biotype=NULL,contrast=NULL,p.value=NULL) standardGeneric("pathways"))




#' @rdname wgcnaDbLite-class
#' @description pathways will return the qusage enrichment pathway descriptions for a given module color.  if you want the entire module enrichment then the Module.color and tx.Biotype must also be a color.  if you want specific traits of interest then set the tx.biotype to the trait.  the contrast is to gather the enrichment data for the comparison either pHSC.vs.LSC (comparison=phsc) or Blast.vs.LSC (comparison=blast)  [lower case only!]
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param tx.Biotype enter hte module color if you want the full module enrichment, or a tx.biotype for a fixed trait
#' @param contrast  either phsc or blast
#' @export
setMethod("pathways", "qusageDbLite", function(x,Module.color=NULL,tx.Biotype=NULL,contrast=NULL,p.value=0.05) {


  ##FIX ME: safety check for the input color using dbListTables,  signature
  allcolors<-dbListTables(dbconn(x))
  stopifnot(Module.color%in%allcolors ==TRUE)

  #drivers are defined as the most significantly cross correlated genes
  if(is.null(contrast)==FALSE){
  sql<-paste0("select pathway_name, logFC, pvalue, FDR, contrastKey from ",Module.color," where colorKey='",Module.color,"'"," and bioKey='",tx.Biotype,"'"," and contrastKey='",contrast,"'")
  }else{
   sql<-paste0("select pathway_name, logFC, pvalue, FDR, contrastKey from ",Module.color," where colorKey='",Module.color,"'"," and bioKey='",tx.Biotype,"'")
 }
  res<-data.frame(dbGetQuery(dbconn(x),sql))
  res<-res[which(res$pathway_name!="NA"),]
  res<-res[which(as.numeric(res$pvalue)<p.value),]
  class(res$pvalue)<-"numeric"
  class(res$logFC)<-"numeric"
  class(res$FDR)<-"numeric"
  return(res)

 })


setGeneric("kexpEnrich",function(x,contrast=NULL,p.value=NULL) standardGeneric("kexpEnrich"))

#' @rdname wgcnaDbLite-class
#' @description kexpEnrich will return the qusage enrichment pathway descriptions for a given entire kexp, not module specific, but every gene in te kexp. The contrast is to gather the enrichment data for the comparison either pHSC.vs.LSC (comparison=phsc) or Blast.vs.LSC (comparison=blast)  [lower case only!]
#' @param x this is the SQLite db
#' @param contrast  either phsc or blast
#' @param p.value  significanse alpha threshold
#' @export
setMethod("kexpEnrich", "qusageDbLite", function(x,contrast=NULL,p.value=0.05) {
 
  if(is.null(contrast)==FALSE){
  sql<-paste0("select * from kexp where colorKey='kexp' and bioKey='kexp' and contrastKey='",contrast,"'")
  }else{
  sql<-paste0("select * from kexp where colorKey='kexp' and bioKey='kexp'")
  }
  res<-as.data.frame(dbGetQuery(dbconn(x),sql))
  res<-res[which(res$p.Value<p.value),]
  colnames(res)<-c("pathway_name","logFC","pvalue","FDR","colorKey","bioKey","contrastKey")
  return(res[,!grepl("row_names",colnames(res))] )

 })


setGeneric("goEnrich",function(x,Module.color=NULL,p.value=NULL) standardGeneric("goEnrich"))

#' @rdname wgcnaDbLite-class
#' @description goEnrich will return the WGCNA experimental enrichment pathway descriptions for a given  gene level kexp, module specific. The contrast is to gather the enrichment data for the comparison is go [lower case only!]
#' @param x this is the SQLite db
#' @param Module.color  this is the gene module to query Goenrich data
#' @param p.value alpha level
#' @export
setMethod("goEnrich", "wgcnaDbLite", function(x,Module.color=NULL,p.value=0.05) {
 sql<-paste0("select * from go where bioKey='kexp' and contrastKey='go' and colorKey='",Module.color,"'")
  res<-as.data.frame(dbGetQuery(dbconn(x),sql))
  res<-res[which(res$p.Value<p.value),]
  pathway.id<-grep("pathway.name",colnames(res))
  p.id<-grep("p.Value",colnames(res))
  f.id<-grep("FDR",colnames(res))
  col.id<-grep("colorKey",colnames(res))
  return(res[,c(pathway.id,p.id,f.id,col.id)] )

 })

setGeneric("listModuleColors",function(x,...) standardGeneric("listModuleColors"))

#' @rdname wgcnaDbLite-class
#' @param x wgcnaDblite SQL-lite database instance
#' @param ... additional parameters releated to annotation database
#' @export
setMethod("listModuleColors", "wgcnaDbLite", function(x, ...) { # {{{
  md <- dbListTables(dbconn(x))
  md<-md[!grepl("metadata",md)]
  md<-md[!grepl("go",md)]
  return(md)
}) # }}}




setGeneric("pickPathway",function(x,p.value=0.05,keyWord=NULL,contrast=NULL) standardGeneric("pickPathway"))




#' @rdname wgcnaDbLite-class
#' @description pickPathway will return the qusage enrichment pathway descriptions for a given module color and a given function that matches the keyword search in the overal module.   the contrast is to gather the enrichment data for the comparison either pHSC.vs.LSC (comparison=phsc) or Blast.vs.LSC (comparison=blast)  [lower case only!]
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param contrast  either phsc or blast
#' @param p.value numeric threshold
#' @param keyWord the name must match what is in the pathways list 
#' @export
setMethod("pickPathway", "qusageDbLite", function(x,p.value=NULL,keyWord=NULL,contrast=NULL) {
  ##FIX ME: safety check for the input color using dbListTables,  signature
  allcolors<-listModuleColors((x))
  df<-list()
  #drivers are defined as the most significantly cross correlated genes
  for(colR in allcolors){  
     if(colR!="kexp"){
    pathd<-pathways((x),Module.color=colR,tx.Biotype=colR,contrast=contrast,p.value=p.value)
    if(nrow(pathd)==0){
    next
    }
    }else if(colR=="kexp"){
  pathd<-kexpEnrich((x),contrast=contrast) 
  pathd<-pathd[,!grepl("colorKey",colnames(pathd))]
  pathd<-pathd[,!grepl("bioKey",colnames(pathd))] 
  }
   ##ranks the absolutre logFC increasing direction
   pathd$ranking<-rank(abs(as.numeric(pathd$logFC)))
   pathd$module.size<-nrow(pathd)
   pathd$query.size<-NA
   pathd$rankTotal<-sum(pathd$ranking)
   pathd$color<-colR
  if(any(grepl(keyWord,pathd$pathway_name,ignore.case=TRUE))){
   ##subset
   res<-pathd[grep(keyWord,pathd$pathway_name,ignore.case=TRUE),]
   marginal<-data.frame(pathway_name=keyWord,logFC=sum(as.numeric(res$logFC)),pvalue=sum(as.numeric(res$pvalue)),FDR=sum(as.numeric(res$FDR)),contrastKey="total",ranking=sum(as.numeric(res$ranking)),module.size=unique(res$module.size),query.size=nrow(res),rankTotal=unique(res$rankTotal),color=colR,row.names="Total")
   query.match<-rbind(res,marginal)
    }else{
# query.match<-data.frame(pathway_name=NA,logFC=NA,pvalue=NA,FDR=NA,contrastKey=NA,ranking=0,module.size=nrow(pathd),query.size=0,rankTotal=unique(pathd$rankTotal))
  next 
    }
  df[[colR]]<-query.match
  }##colR loop
   
  return(df)

 })


###FIX ME: add a consensus pathways function
