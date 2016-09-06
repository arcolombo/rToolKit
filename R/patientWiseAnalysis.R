#' @title patient Differential Expression analysis
#' @description patientPlot* is exploratory methods for expression patterns, this method is downstream of a selected DE caller. the inputs are a repeat stage level kexp, and the DE object from either cqnDE, poissonFit, edgeR, or limma calls.  The topSelected repeats are used to plot the overall heatmap, and then heatmap each tx_biotype separately.   the topSelected are then categorized based on ordinal groups, and the group frequencies are plotted.  pvalue histograms are included.  delta plots are included (logFC). 
#' @import pvclust
#' @import dendsort
#' @import arkas
#' @import ComplexHeatmap
#' @param kexp for DE poissonSeq kexp should repeat stage, for cqn kexp should be stage transcrpit
#' @param whichDE the de caller cqn, poissonSeq limma , limma should call stageWiseaNALYSIS
#' @param comparison   character the comparison used for 2 group cmpr/control
#' @param control  the character in the denominator of 2 group DE call
#' @param whichDelta  characther  either LSC/pHSC or Blast/LSC
#' @param patientID   character  ID to subset
#' @param orderDE    boolean ,  if true the TPMs of the topNamed repeat tags are plotted,  if FALSE the TPMs are all plotted for global view
#' @export 
#' @return a list of all the monotonic ordinal groups  
patientWiseAnalysis<-function(kexp,whichDE=c("cqn","poissonSeq","limma"),comparison="LSC",control="pHSC",whichDelta=c("delta1","delta2"),patientID="SU353", orderDE=TRUE){

whichDE<-match.arg(whichDE,c("cqn","poissonSeq","limma"))
  whichDelta<-match.arg(whichDelta,c("delta1","delta2"))
  if(whichDE=="poissonSeq"){
##poissonfit  colnames(tt) = c("names","logfc", fdr)
  pf<-poissonFit(kexp,comparison=comparison,control=control,whichDelta=whichDelta)
   patient<-kexpByPatient(kexp,patientID=patientID)
   ##global heat of all poissonSeq selected
   message("heatmap for patient: ",patientID)
   drawDF(patient,res=pf, fromList=TRUE) 
   message("DE composition for patient: ",patientID)
   #plotFrequency(patient,topDE=pf,whichDelta="delta1",isAdjusted=TRUE)
   plot.new()
   hist(pf$fdr,main=paste0("FDR values for patient: ",patientID),xlab="FDR")    
 ## ordinal grouping of topNames @pvalue
   ####call plotMonotonicity
   stageId1<-grep("pHSC",colnames(patient))
   stageId2<-grep("LSC",colnames(patient))
   stageId3<-grep("Blast",colnames(patient))
   mt<-collapseTpm(patient,"tx_id") 
   if(orderDE==TRUE){
   message("subsetting repeat tpm ",patientID," by DE selection")
   mt<-mt[rownames(mt)%in%pf$names,]
   }
    message("pHSC>Blast>LSC at pvalue: 0.05")
   
   pHSC.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId2,globalMid=stageId3)
   plotMonotonicity(patient,
                    nominalGroup=pHSC.max.min.LSC,
                    stageId1=stageId1,
                    stageId2=stageId2,
                    stageId3=stageId3,
                    mt, 
                    read.cutoff=1,
                    stage3="Blast", 
                    Max=colnames(patient)[stageId1],
                    Mid=colnames(patient)[stageId3],
                    Min=colnames(patient)[stageId2] )
   ## patient delta plots
  blast.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId3,
                     globalMin=stageId2,
                     globalMid=stageId1)
   plotMonotonicity(patient,
                    nominalGroup=blast.max.min.LSC,
                    stageId1=stageId1,
                    stageId2=stageId2,
                    stageId3=stageId3,
                    mt,
                    read.cutoff=1,
                    stage3=stage3, 
                    Max=colnames(patient)[stageId3],
                    Mid=colnames(patient)[stageId1],
                    Min=colnames(patient)[stageId2])

 LSC.max.min.pHSC<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId1,globalMid=stageId3)
 plotMonotonicity(patient,
                 nominalGroup=LSC.max.min.pHSC,
                 stageId1=stageId1,
                 stageId2=stageId2,
                 stageId3=stageId3,
                 mt,
                 read.cutoff=1,
                 stage3=stage3,
                 Max=colnames(patient)[stageId2],
                 Mid=colnames(patient)[stageId3], 
                 Min=colnames(patient)[stageId1] )

  ### 2>1>3
  LSC.max.min.blast<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId3,globalMid=stageId1)
 plotMonotonicity(patient,
                  nominalGroup=LSC.max.min.blast,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1,
                  stage3=stage3, 
                  Max=colnames(patient)[stageId2], 
                  Mid=colnames(patient)[stageId1],
                  Min=colnames(patient)[stageId3])


   pHSC.to.blast<-orderMonotonicity(mt,globalMax=stageId3,globalMin=stageId1,globalMid=stageId2)
  if(nrow(pHSC.to.blast)>3){
  plotMonotonicity(patient,
                  nominalGroup=pHSC.to.blast,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1,
                  stage3=stage3, 
                  Max=colnames(patient)[stageId3],
                  Mid=colnames(patient)[stageId2],
                  Min=colnames(patient)[stageId1])
   } else {
    print(pHSC.to.blast)
   }
  ## 1>2>3
  
  blast.to.pHSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId3,globalMid=stageId2)
  if(nrow(blast.to.pHSC)>3) {
  plotMonotonicity(patient,
                  nominalGroup=blast.to.pHSC,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1, 
                  stage3=stage3, 
                  Max=colnames(patient)[stageId1],
                  Mid=colnames(patient)[stageId2],
                  Min=colnames(patient)[stageId3])

   } else {
   print(blast.to.pHSC)
   }



  } else if(whichDE=="cqn"){
### cqn list names = "fitted","topTags", "RPKM.cqn.log2"
  message("for cqn option please use the full stage kexp")
  cqn.full<-cqnDE(fullKexp,comparison=comparison,control=control)
  Repeat.counts<-cqn.full[["RPKM.cqn.log2"]]
  Repeat.counts<-Repeat.counts^2
  Repeat.counts<-Repeat.counts[!grepl("^ERCC",rownames(Repeat.counts)),]
  Repeat.counts<-Repeat.counts[!grepl("^ENST",rownames(Repeat.counts)),]
  fullTags<-rownames(cqn.full["topTags"][[1]])
  fullDE<-as.data.frame(cqn.full["topTags"][1])
  colnames(fullDE)<-c("length","gccontent","logFC","logCPM","LR","PValue","FDR")
  res<-data.frame(names=rownames(cqn.full$topTags),
                  fdr=fullDE$FDR ,
                  stringsAsFactors=FALSE)
   patient<-kexpByPatient(findRepeats(kexp),patientID=patientID)
   ##global heat of all poissonSeq selected
   message("heatmap for patient: ",patientID)
   drawDF(patient,res=res, fromList=TRUE)
   message("DE composition for patient: ",patientID)
   #plotFrequency(patient,topDE=pf,whichDelta="delta1",isAdjusted=TRUE)
   plot.new()
   hist(res$fdr,main=paste0("FDR values for patient: ",patientID),xlab="FDR")
 ## ordinal grouping of topNames @pvalue
   ####call plotMonotonicity
   stageId1<-grep("pHSC",colnames(patient))
   stageId2<-grep("LSC",colnames(patient))
   stageId3<-grep("Blast",colnames(patient))
   mt<-collapseTpm(patient,"tx_id")
   if(orderDE==TRUE){
   message("subsetting repeat tpm ",patientID," by DE selection")
   mt<-mt[rownames(mt)%in%res$names,]
   }
    message("pHSC>Blast>LSC at pvalue: 0.05")
   pHSC.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId2,globalMid=stageId3)
   plotMonotonicity(patient,
                    nominalGroup=pHSC.max.min.LSC,
                    stageId1=stageId1,
                    stageId2=stageId2,
                    stageId3=stageId3,
                    mt,
                    read.cutoff=1,
                    stage3="Blast",
                    Max=colnames(patient)[stageId1],
                    Mid=colnames(patient)[stageId3],
                    Min=colnames(patient)[stageId2] )
   ## patient delta plots
  message("Blast>pHSC>LSC at pvalue 0.05")
  blast.max.min.LSC<-orderMonotonicity(mt,globalMax=stageId3,
                     globalMin=stageId2,
                     globalMid=stageId1)
   plotMonotonicity(patient,
                    nominalGroup=blast.max.min.LSC,
                    stageId1=stageId1,
                    stageId2=stageId2,
                    stageId3=stageId3,
                    mt,
                    read.cutoff=1,
                    stage3=stage3,
                    Max=colnames(patient)[stageId3],
                    Mid=colnames(patient)[stageId1],
                    Min=colnames(patient)[stageId2])

 message("LSC > Blast>pHSC pvalue 0.05")
 LSC.max.min.pHSC<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId1,globalMid=stageId3)
 plotMonotonicity(patient,
                 nominalGroup=LSC.max.min.pHSC,
                 stageId1=stageId1,
                 stageId2=stageId2,
                 stageId3=stageId3,
                 mt,
                 read.cutoff=1,
                 stage3=stage3,
                 Max=colnames(patient)[stageId2],
                 Mid=colnames(patient)[stageId3],
                 Min=colnames(patient)[stageId1] )


 ### 2>1>3
  LSC.max.min.blast<-orderMonotonicity(mt,globalMax=stageId2,globalMin=stageId3,globalMid=stageId1)
 plotMonotonicity(patient,
                  nominalGroup=LSC.max.min.blast,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1,
                  stage3=stage3,
                  Max=colnames(patient)[stageId2],
                  Mid=colnames(patient)[stageId1],
                  Min=colnames(patient)[stageId3])


   pHSC.to.blast<-orderMonotonicity(mt,globalMax=stageId3,globalMin=stageId1,globalMid=stageId2)
  if(nrow(pHSC.to.blast)>3){
  plotMonotonicity(patient,
                  nominalGroup=pHSC.to.blast,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1,
                  stage3=stage3,
                  Max=colnames(patient)[stageId3],
                  Mid=colnames(patient)[stageId2],
                  Min=colnames(patient)[stageId1])
 
  } else {
    print(pHSC.to.blast)
   }

 blast.to.pHSC<-orderMonotonicity(mt,globalMax=stageId1,globalMin=stageId3,globalMid=stageId2)
  if(nrow(blast.to.pHSC)>3) {
  plotMonotonicity(patient,
                  nominalGroup=blast.to.pHSC,
                  stageId1=stageId1,
                  stageId2=stageId2,
                  stageId3=stageId3,
                  mt,
                  read.cutoff=1,
                  stage3=stage3,
                  Max=colnames(patient)[stageId1],
                  Mid=colnames(patient)[stageId2],
                  Min=colnames(patient)[stageId3])

   } else {
   print(blast.to.pHSC)
   }



  } else {
  ###limma call from arkas ruvRWA


  }
 
 message(paste0("the ordinal group elements for ",nrow(mt)," ",
         nrow(pHSC.max.min.LSC),"+",
         nrow(blast.max.min.LSC),"+",
         nrow(LSC.max.min.pHSC),"+",
         nrow(pHSC.to.blast),"+",
         nrow(blast.to.pHSC), " =",
         sum(nrow(pHSC.max.min.LSC),nrow(blast.max.min.LSC),nrow(LSC.max.min.pHSC),nrow(pHSC.to.blast), nrow(blast.to.pHSC))))



##FIX ME add other groups  A=B>C  A>B=C  B=C>A 



 return(list(A= pHSC.max.min.LSC,
              B= blast.max.min.LSC,
              C= LSC.max.min.pHSC,
              D= pHSC.to.blast,
              E =  blast.to.pHSC))
 
} ##{{main

