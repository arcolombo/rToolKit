#' selects the top MAD, this was authored by Tim Triche Jr.
#' @param x a matrix
#' @param k the number of top MADs to select
#' @importFrom matrixStats rowMads
#' @export
byMad <- function(x, k=500) { 
  mads <- matrixStats::rowMads(x, na.rm=TRUE)
  x[rev(order(mads))[seq_len(k)],]
}

