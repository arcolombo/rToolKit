#' plots a heatmap of patient specific information using a heatmap and beeswarm for the scatter.  each plot shows fixed patient and fixed repeat biotype
#' @param kexp a kallisto Experiment
#' @param patientID a patient id to split a kexp by
#' @param repeatType a repeat biotype of interest, 'LTR', 'SINE', 'element', etc
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom graphics grid
#' @import beeswarm
#' @importFrom EDASeq plotPCA plotRLE
#' @export
#' @return returns NULL but plots cool stuff
patientPlot<-function(kexp,patientID=NULL,repeatType=NULL){

 readkey<-function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

#subset repeats
repKexp<-findRepeats(kexp)
##grab patient

patient.kexp<-kexpByPatient(repKexp,patientID=patientID)

##grab repeat gene biotype
if(is.null(repeatType)==FALSE){
patient.type<-kexpByType(patient.kexp,biotype=repeatType)
  #add complex heatmap  plots
 selectK<-nrow(patient.type)/2
 typed<-Heatmap(byMad(counts(patient.type),k=selectK),column_title=paste0(patientID," ",repeatType," byMad Top Half LTR CPM"))
  draw(typed)
  readkey()
   } else {
 patient.ltr<-kexpByType(patient.kexp,biotype="LTR")
 patient.sine<-kexpByType(patient.kexp,biotype="SINE")
 #add complex heatmap  plots  ###
 ##FIX ME :: add mutation information as a annotation
 ##add annotation heatmaps
  selectK<-nrow(patient.ltr)/2
  ltr<-Heatmap(byMad(counts(patient.ltr),k=selectK),column_title=paste0(patientID," byMad Top Half LTR CPM"))
  draw(ltr) 
  readkey()
  selectM<-nrow(patient.sine)
  sine<-Heatmap(byMad(counts(patient.sine),k=selectM),column_title=paste0(patientID," byMad All SINE CPM"))
   draw(sine)
   readkey()
  ##bySd
 ltr.sd<-Heatmap(bySd(counts(patient.ltr),k=nrow(patient.ltr)/2),column_title=paste0(patientID," bySd Top Half LTR CPM"))
  draw(ltr.sd) 
  readkey()
  sine.sd<-Heatmap(bySd(counts(patient.sine),k=nrow(patient.sine)),column_title=paste0(patientID," bySd All SINE CPM"))
   draw(sine.sd)
   readkey()
  ##add beeswarm plots

  #plots top 20 byMAD individual repeats across clonal evolution
  beePatient(patient.ltr,patientID=patientID,repeat_biotype="LTR",selected=20)
  beePatient(patient.sine,patientID=patientID,repeat_biotype="SINE",selected=20)
  
  #now plot all repeats together
  plotRLE(counts(patient.ltr))
  title("All LTR CPM")
  readkey()
  plotRLE(counts(patient.sine))
  title("All SINEs CPM")
  readkey()
  
  #plot by columns
  beeColumns(patient.ltr,patientID=patientID,repeat_biotype="LTR",selected=nrow(patient.ltr),what="count")
  readkey() 
  beeColumns(patient.sine,patientID=patientID,repeat_biotype="SINE",selected=nrow(patient.sine),what="count") 
  readkey()
  beeColumns(patient.ltr,patientID=patientID,repeat_biotype="LTR",selected=nrow(patient.ltr),what="tpm") 
  readkey()
  beeColumns(patient.sine,patientID=patientID,repeat_biotype="SINE",selected=nrow(patient.sine),what="tpm")
 readkey()
  ##calculate pearson coef   cor() 
  #FIX ME: add mutation annotations plotRLE(counts(SU583.ltr))
   ##FIX ME::  calculate TMM and normalize by library size

   ##FIX ME ::  run CQN   



}




} #{{{ main
