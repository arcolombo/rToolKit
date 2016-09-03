#' @title runs RUV to normalize stage wise expression 
#' @description plots the heatmap of top byMAD repeats. shows the composition of biotypes of the significant repeats, plots the delta beeswarm plots of each stage to show a global trend of the changes among states.
#' @param stageKexp kexp at the stages
#' @param p.value numeric significant cutoff
#' @importFrom RUVSeq RUVg
#' @import ComplexHeatmap
#' @import grid
ruvStageWiseAnalysis<-function(stageKexp,p.value=0.05,cutoff=1){


  design<-model.matrix(~0+sapply(strsplit(pData(stageKexp)$ID,"_"),function(x) x[1]))
  metadata(stageKexp)$design<-design
  weights<-ruvNormalization(stageKexp,byLevel="tx_id")


} #{{{ main
