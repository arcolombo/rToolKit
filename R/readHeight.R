#' read key pause for plotting
#' @export
 readHeight<-function(){
    cat("Enter height to cut:\n")
    line <- readline()
    message(paste0("Cutting :",line))
    return(as.numeric(line))
}


