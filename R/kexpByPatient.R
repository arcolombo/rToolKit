#' splits kexp by patient ID
#' @param kexp a kallisto experiment
#' @param patientID a patient id to split columns
#' @export 
kexpByPatient <-function(kexp,patientID=NULL){

return(kexp[,grep(patientID,colnames(kexp),ignore.case=TRUE)])
}
