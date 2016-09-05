#' @title omits SU583 from kexp
#' @description SU583 donor was not taken the same as the others but was taken post relapse, this sample should be omitted.
#' @param kexp a kexp at any level with SU583 as column
#' @export
#' @return a kexp
omitSU583<-function(kexp){
return(kexp[,!grepl("SU583",colnames(kexp))])
}
