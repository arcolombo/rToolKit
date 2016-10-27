#' @title omits SU583 from kexp
#' @description SU583 donor was not taken the same as the others but was taken post relapse, this sample should be omitted.
#' @param kexp a kexp at any level with SU583 as column
#' @param relapsedID  character SU351 and SU583 are reported relapsed, removing
#' @export
#' @return a kexp
omitRelapsed<-function(kexp,relapsedID=c("SU351","SU583")){
message("ommitting SU351 and SU583 due to relapse condition")
kexp<-kexp[,!grepl("SU583",colnames(kexp))]
kexp<-kexp[,!grepl("SU351",colnames(kexp))]

return(kexp)
}
