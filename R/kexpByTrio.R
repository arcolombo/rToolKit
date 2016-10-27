#' @title this will subset a kexp by patient trio data
#' @description given a trio time series we subset a kexp so that each patient trio exists across all time points
#' @export
kexpByTrio<-function(kexp,stage1="pHSC",stage2="LSC",stage3="Blast"){

  p.id<-sapply(strsplit(colnames(kexp),"_"),function(x) x[2])
 trio.flag<-sapply(p.id,function(x) which( length(grep(x,colnames(kexp)))==3) )
  pass.flag<-sapply(trio.flag,function(x) length(x>0))
  names.pass<-unique(names(pass.flag[which(pass.flag>0)]))
  t<-c(paste0(stage1,"_",names.pass),paste0(stage2,"_",names.pass),paste0(stage3,"_",names.pass))
  return(kexp[,colnames(kexp)%in%t])
 }
