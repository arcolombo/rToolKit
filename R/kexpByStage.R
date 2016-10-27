#' splits a kallisto experiment by stages
#' @param kexp a kallisto Experiment
#' @param stage1 first stage
#' @param stage2 second stage
#' @param stage3 third stage
#' @export
kexpByStage<-function(kexp,stage1="pHSC",stage2="LSC",stage3="Blast"){

 id1<-grep(stage1,colnames(kexp))
 id2<-grep(stage2,colnames(kexp))
 id3<-grep(stage3,colnames(kexp))
 
return(kexp[,c(id1,id2,id3)])
}
