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



setGeneric("driversBy",function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL,...) standardGeneric("driversBy"))



#' @rdname wgcnaDbLite-class
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param trait  this is the trait of interest 
#' @param p.value p.value significance threshold
#' @export
setMethod("driversBy", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,p.value=NULL) {
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
  gcor<-paste0("GCor_",trait)
  pval.Trait<-paste0(p.prefix,trait)
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

setGeneric("traitsBy",function(x,Module.color=NULL,p.value.type="studentT",p.value=NULL,...) standardGeneric("traitsBy"))


#' @rdname wgcnaDbLite-class
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param p.value p.value significance threshold
#' @export
setMethod("traitsBy", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",p.value=0.05) {
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
  
  #trait was not specified so return all traits 
  sql<-paste0("select * from ",Module.color," where colorKey='",Module.color,"'")

   res<-as.data.frame(dbGetQuery(dbconn(x),sql),
                    stringsAsFactors=FALSE)
   ##only return one specified p.value type 
   id<-!grepl(paste0("^",alt.p.prefix,"*"),colnames(res))
   #select one module and every trait (column)
   res<-res[,id]
 
    
   ###  loop through each trait and find drivers per trait
   outList<-list()
   traits.id<-grep("^GCor_",colnames(res))
   module.id<-grep(paste0("MM",Module.color),colnames(res))
   rownames.id<-grep("row_names",colnames(res))
   for(i in 1:length(traits.id)){
   col.id<-which(substring(colnames(res)[traits.id[i]],6)==substring(colnames(res),8))
    out<-res[,c(rownames.id,module.id,traits.id[i],col.id)] 
   p.col<-paste0(p.prefix,substring(colnames(res)[traits.id[i]],6))
    out<-out[which(out[,grepl(p.col,colnames(out))]<=p.value),]
    if(nrow(out)>0){
    #only store non-empty sets
    outList[[i]]<-out
    names(outList)[i]<-paste0(Module.color,"_",substring(colnames(res)[traits.id[i]],6))
    } #if
  }#for
  
    return(outList)

 })



### FIX ME: add verbose plots to plot significant correlations consensus




