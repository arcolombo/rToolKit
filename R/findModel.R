#' @title find the repeats falling into orders
#' @description  find the repeats fallingo into clonal ordering with the alternative model and hypothesis model input for element selection. this has no p value association but expression pattern overview 
#' @param kexp a kallistoExperiment of a clonal patient kexp repeat
#' @param hypothesis  an ordering character
#' @param alternative an alternative ordering character
#' @param patientID  a patient character identifier
#' @param repeatType a repeat character biotype
#' @param read.cutoff integer a min cutoff 
#' @param outputReport boolean if true print a CSV
#' @param outputDir   character path
#' @importFrom graphics pie
#' @export
#' @return a repeat listing for each hypothesis 
findModel<-function(kexp,stage1="pHSC", stage2="LSC",stage3="Blast",patientID=NULL,repeatType=NULL,read.cutoff=1,outputReport=FALSE,outputDir=NULL, what=c("counts","tpm")){
  what<-match.arg(what,c("counts","tpm"))
 #the input is a clonal kexp, findModel will split
  if(is.null(patientID)==FALSE){
 kexp<-kexpByPatient(kexp,patientID=patientID)
  }
  stageId1<-grep(stage1,colnames(kexp))
  stageId2<-grep(stage2,colnames(kexp))
  stageId3<-grep(stage3,colnames(kexp))

  if(what=="counts"){
  mt<-collapseBundles(kexp,"tx_id",read.cutoff=read.cutoff)
  } else {
  mt<-collapseTpm(kexp,"tx_id",read.cutoff=read.cutoff)
  }
  ##FIX ME : add option for normalization
  
  ##fix me add all query groups
  ## Non-Monotonic Group Strict Ordering only
  ## 1>3>2 
  pHSC.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId2,globalMid=stageId3)
  plotMonotonicity(kexp,nominalGroup=pHSC.max.min.LSC,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3,Max=colnames(kexp)[stageId1],Mid=colnames(kexp)[stageId3],Min=colnames(kexp)[stageId2] )
  ### 3>1>2
   blast.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId3,globalMin=stageId2,globalMid=stageId1)
   plotMonotonicity(kexp,nominalGroup=blast.max.min.LSC,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3, Max=colnames(kexp)[stageId3],Mid=colnames(kexp)[stageId1],Min=colnames(kexp)[stageId2])

  ###2>3>1
  LSC.max.min.pHSC<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId1,globalMid=stageId3)
 plotMonotonicity(kexp,nominalGroup=LSC.max.min.pHSC,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3,Max=colnames(kexp)[stageId2],Mid=colnames(kexp)[stageId3], Min=colnames(kexp)[stageId1] )

  ### 2>1>3
  LSC.max.min.blast<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId3,globalMid=stageId1)
 plotMonotonicity(kexp,nominalGroup=LSC.max.min.blast,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3, Max=colnames(kexp)[stageId2], Mid=colnames(kexp)[stageId1], Min=colnames(kexp)[stageId3])

  ###
  
  ### Mono Tonic Group
  ## 1<2<3
  pHSC.to.blast<-orderMonotonicity(mt,globalMax=stageId3,globalMin=stageId1,globalMid=stageId2)
 plotMonotonicity(kexp,nominalGroup=pHSC.to.blast,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3, Max=colnames(kexp)[stageId3],Mid=colnames(kexp)[stageId2],Min=colnames(kexp)[stageId1])

  ## 1>2>3
  blast.to.pHSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId3,globalMid=stageId2)
 plotMonotonicity(kexp,nominalGroup=blast.to.pHSC,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=read.cutoff,stage3=stage3, Max=colnames(kexp)[stageId1],Mid=colnames(kexp)[stageId2],Min=colnames(kexp)[stageId3])

  ####

  ####Non Decreasing Group
  ##need to specif
 
  
   
   ##FIX ME:  print out the top 40-50 repeats
 ##print out

 if(outputReport==TRUE){
 write.csv(df2,file=paste0(outDir,"/",stage1,"_",stage3,"_",stage2,".csv",row.names=TRUE))
  ##FIX ME:  create a pdf with every repeat 1 pdf
  ##save all data
 } else {
  return(NULL) 
 }



  ##FIX ME:  should not the groups  A > B > C   and A < B < C  be equal to N? 
  ###FIX ME ::  N(pHSC>LSC) + N(LSC>pHSC ) = N_total
} #{{{main
