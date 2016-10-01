#' @rdname wgcnaDbLite-class
#' @export
setMethod("dbconn", "wgcnaDbLite", function(x) return(x@con))



#' @rdname wgcnaDbLite-class
#' @param object wgcnaDbLite SQL-lite annotation database
#' @import ensembldb
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
setMethod("metadata", "wgcnaDbLite", function(x, ...) { # {{{
  md <- dbGetQuery(dbconn(x), "select * from metadata")
  rownames(md) <- md$name
  return(md)
}) # }}}


#FIX ME: add correlation queries of the database
###for colorname  select from table color where colorKey=color 
###  


###FIX ME:  add filter scripts where one can filter based on p.value


### FIX ME: add verbose plots to plot significant correlations as potential drivers

### FIX ME: add Cormap functions .




