#' Constructor for an qusageDbLite object.
#'
#' @param x the name of the sqlite file
#' @param ... any additional arguments relating to ENSEMBL annotation database.
#' @return an qusageDbLite object.
#' 
#' @export
qusageDbLite <- function (x, ...) {
  qusagedb <- wgcnaDbLite(x, ...)
  class(qusagedb) <- "qusageDbLite"
  return(qusagedb)
}

