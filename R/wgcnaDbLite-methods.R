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


setGeneric("modulesBy",function(x,Module.color=NULL,p.value.type="studentT",trait=NULL,...) standardGeneric("modulesBy")) 
 
  
#' @rdname wgcnaDbLite-class
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param trait  this is the trait of interest 
#' @export
#' @import TxDbLite
setMethod("modulesBy", "TxDbLite", function(x,Module.color=NULL,p.value.type="studentT",trait=NULL) {
  p.value.type<-match.arg(p.value.type,c("studentT","fisher"))
  if(p.value.type=="studentT"){
  p.prefix<-"p_GCor_"
  alt.p.prefix<-"pf_GCor_"
  } else{
  p.prefix<-"pf_GCr_"
  alt.p.prefix<-"p_GCor_"
  }
  ## safety check for the input color using dbListTables
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
   }
    return(res)
  
 }) 


#' @rdname wgcnaDbLite-class
#' @param x this is the SQLite db
#' @param Module.color this is a selected table in the db, a module color 
#' @param p.value.type either studentT or fisher, two types of p.value tests were done.
#' @param trait  this is the trait of interest 
#' @export
#' @import TxDbLite
setMethod("modulesBy", "wgcnaDbLite", function(x,Module.color=NULL,p.value.type="studentT",trait=NULL) {
  p.value.type<-match.arg(p.value.type,c("studentT","fisher"))
  if(p.value.type=="studentT"){
  p.prefix<-"p_GCor_"
  alt.p.prefix<-"pf_GCor_"
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
   }
    return(res)
  
 })


###FIX ME:  add filter scripts where one can filter based on p.value


### FIX ME: add verbose plots to plot significant correlations as potential drivers

### FIX ME: add Cormap functions .




