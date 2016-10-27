#' Create a wgcnaDbLite object (usually called by subclass constructors)
#' 
#' @param x an annotation instance of wgcnaDbLite 
#' @param path  where it lives (default: ".")
#' @param ...   additional parameters 
#' @return   a TxDbLite object 
#' @importFrom DBI dbDriver dbConnect dbListTables dbGetQuery
#' @export
qusageDbLite <- function(x, path=".", ...) {
  options(useFancyQuotes = FALSE)
  if (!grepl(".sqlite$", x)) x <- paste0(x, ".sqlite")
  x <- paste(path, x, sep=.Platform$file.sep)
  lite <- dbDriver("SQLite")
  con <- dbConnect(lite, dbname = x, flags = SQLITE_RO)
  tables <- dbListTables(con)
  names(tables) <- tables
  Tables <- lapply(tables, function(x)
    colnames(dbGetQuery(con, paste0("select * from ", x, " limit 1"))))
  new("qusageDbLite", con = con, tables = Tables)
}

