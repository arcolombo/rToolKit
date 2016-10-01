#' A parent class for wgcna databases including annotations
#' @rdname wgcnaDbLite-class
#' @import methods
#' @slot con a DBI connection
#' @slot tables a list of annotation entries
#' @export 
setClass("wgcnaDbLite",
  representation(con="DBIConnection",tables="list"),
  prototype=list(con=NULL,tables=list())
)

#' A parent class for Ensembl based network annotations
#' @rdname wgcnaDbLite-class
#' @export 
setClass("wgcnaDbLite",contains="TxDbLite")

