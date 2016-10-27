#' @title plots the monotonic order
#' @description plots the monotonic ordering first users must call orderMonotonicty and the nominalGroup is the results from orderMontonicity
#' @param kexp a kallisto Experiment patient repeat level
#' @param nominalGroup a monotonic order
#' @param stageId1  the first clone stage
#' @param stageId2  the second clone stage
#' @param stageId3  the third stage
#' @param mt  a collapse bundles matrix
#' @param read.cutoff the cutoff for CPM min floor
#' @param Max a column of a kexp
#' @param Mid a column of a kexp
#' @param Min a column of a kexp
#' @import arkas
#' @import graphics
#' @import beeswarm
#' @export
#' @return some plots 
plotMonotonicity<-function(kexp,nominalGroup=NULL,stageId1=stageId1,stageId2=stageId2,stageId3=stageId3,mt,read.cutoff=1, stage3=stage3,Max=globalMax,Mid=globalMid,Min=globalMin){

hypothesis<-nominalGroup
  ##find annotations
  df<-data.frame(names=rownames(hypothesis),
                 tx_biotype=rowRanges(kexp)[rownames(hypothesis)]$tx_biotype,
                 gene_biotype=rowRanges(kexp)[rownames(hypothesis)]$gene_biotype)

  df2<-cbind(hypothesis,df)

  dTab<-(table(df[,2]))
  #####freq bar chart of nominal category variables
  barplot(dTab,las=2,cex.names=0.6,main=paste0("TxDb Freq N=",nrow(hypothesis),"/",nrow(mt)," ",Max,">",Mid,">",Min  )  )
  
  readkey()

  ######## relative freq bar chart of category variables
   barplot((dTab/sum(dTab))*100,las=2,cex.names=0.6, main=paste0("TxDb Freq N=",nrow(hypothesis),"/",nrow(mt)," ", Max,">",Mid,">",Min  ))
  readkey()


  ######gene biotype
  barplot(table(df[,3]),las=2,cex.names=0.6, main=paste0("TxDb Freq ",Max,">",Mid,">",Min  ) )


  #########pie chart 
  pie(dTab[which(dTab>=read.cutoff)] )
  dTabFreq <- prop.table(dTab[which(dTab>=read.cutoff)] )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("TxBiotypes N=",nrow(hypothesis),"/",nrow(mt)," ",Max,">",Mid,">",Min,">",read.cutoff ))
  readkey()


 dTab<-(table(df[,3]))
  pie(dTab)
  dTabFreq <- prop.table(dTab )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("GeneBiotypes N=",nrow(hypothesis),"/",nrow(mt)," ",Max,">",Mid,">",Min,">0"))
  readkey()


  boxplot(asinh(df2[,c(stageId1,stageId2,stageId3)]),medcol="red",boxcol="red",whiskcol="red")
   axis(side=4)
  par(new=TRUE)
  beeswarm(asinh(df2[,c(stageId1,stageId2,stageId3)]), pch=16,xlab=paste0(Max,">",Mid,">",Min,">",read.cutoff) )
    readkey()


} #{{{ main
