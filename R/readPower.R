#' read key pause for plotting powers
#' @export
 readPower<-function(){
    cat("Enter desired power :\n")
    line <- readline()
    message(paste0("Power Selected :",line))
    return(as.numeric(line))
}


