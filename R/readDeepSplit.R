#' read key pause for plotting powers
#' @export
 readDeepSplit<-function(){
    cat("Enter desired split level :\n")
    line <- readline()
    message(paste0("Deep Split Selected :",line))
    return(as.numeric(line))
}


