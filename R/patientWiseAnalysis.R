#' @title patient Differential Expression analysis
#' @description patientPlot* is exploratory methods for expression patterns, this method is downstream of a selected DE caller. the inputs are a repeat stage level kexp, and the DE object from either cqnDE, poissonFit, edgeR, or limma calls.  The topSelected repeats are used to plot the overall heatmap, and then heatmap each tx_biotype separately.   the topSelected are then categorized based on ordinal groups, and the group frequencies are plotted.  pvalue histograms are included.  delta plots are included (logFC). 
#' @import pvclust
#' @import dendsort
#' @import arkas
#' @import ComplexHeatmap
patientWiseAnalysis<-function(kexp,topNames=NULL,topDE=NULL,whichDE=c("cqn","poissonSeq","limma"),comparison="LSC",control="pHSC",whichDelta=c("delta1","delta2"),patientID="SU353"){
stopifnot(is.null(topDE)==FALSE) ##must have DE called
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
   plotFrequency(patient,topDE=pf,whichDelta="delta1",isAdjusted=TRUE)
   plot.new()
   hist(pf$fdr,main=paste0("FDR values for patient: ",patientID),xlab="FDR")    
 ## ordinal grouping of topNames @pvalue
    
   ##  delta plots

  } else if(whichDE=="cqn"){
### cqn list names = "fitted","topTags", "RPKM.cqn.log2"

  } else {
  ###limma call from arkas ruvRWA


  }
} ##{{main

