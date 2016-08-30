#' @title plot frequency of biotypes
#' @description plots the histogram and bar charts of biotypes
#' @param kexp  a kallistoExperiment
#' @param topNames names of the top DE list of PVal or adj
#' @param stage1 first stage
#' @param stage2 second stage
#' @param stage3 third stage
#' @param topDE differential expr list
#' @param whichDelta significant DE is a change between states, so which state is being measured?
#' @param p.cutoff numeric p value threshold
#' @param isAdjusted boolean if true than the adj.P.Values are used and labels with show accordingly
#' @export
#' @return cool images
plotFrequency<-function(kexp,topNames=NULL,stage1="pHSC",stage2="LSC",stage3="Blast",topDE=NULL,whichDelta=c("delta1","delta2"),p.cutoff=0.05,isAdjusted=FALSE){

 whichDelta<-match.arg(whichDelta,c("delta1","delta2"))
df<-data.frame(names=rownames(topDE),
               tx_biotype=rowRanges(kexp)[rownames(topDE)]$tx_biotype,
               gene_biotype=rowRanges(kexp)[rownames(topDE)]$gene_biotype)

  dTab<-table(df$tx_biotype)
 if(whichDelta=="delta1") {
  #Delta 1   LSC-pHSC
  if(isAdjusted==FALSE){ 
 barplot(dTab,las=2,cex.names=0.6,main=expression(paste(Delta,"(pHSC,LSC) Top DE Tx Biotype")),ylab=paste0(" Freq N=",nrow(topDE),"/",nrow(kexp)  ),xlab=paste0("p.value=",p.cutoff)  )
  readkey()
  } else {
 barplot(dTab,las=2,cex.names=0.6,main=expression(paste(Delta,"(pHSC,LSC) Top DE Tx Biotype BH")),ylab=paste0(" Freq N=",nrow(topDE),"/",nrow(kexp)  ),xlab=paste0("adj.p.value=",p.cutoff)  )
 }

 #########pie chart 
  pie(dTab )
  dTabFreq <- prop.table(dTab )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("DE TxBiotypes N=",nrow(topDE),"/",nrow(kexp)," ",p.cutoff ))
  readkey()


 dTab<-(table(df$gene_biotype))
  pie(dTab)
  dTabFreq <- prop.table(dTab )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("DE GeneBiotypes N=",nrow(topDE),"/",nrow(kexp)," ",p.cutoff))
  readkey()

  

  } else {
   #for Blast-LSC Delta2
  if(isAdjusted==FALSE){
  barplot(dTab,las=2,cex.names=0.6,main=expression(paste(Delta,"(LSC,Blast) Top DE Tx Biotype")),ylab=paste0(" Freq N=",nrow(topDE),"/",nrow(kexp)),xlab=paste0("p.value=",p.cutoff)  )
  readkey()
  } else {
  barplot(dTab,las=2,cex.names=0.6,main=expression(paste(Delta,"(LSC,Blast) Top DE Tx BiotypeBH")),ylab=paste0(" Freq N=",nrow(topDE),"/",nrow(kexp)),xlab=paste0("adj.p.value=",p.cutoff)  )
   }
  pie(dTab )
  dTabFreq <- prop.table(dTab )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("DE TxBiotypes N=",nrow(topDE),"/",nrow(kexp)," ",p.cutoff ))
  readkey()


 dTab<-(table(df$gene_biotype))
  pie(dTab)
  dTabFreq <- prop.table(dTab )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("DE GeneBiotypes N=",nrow(topDE),"/",nrow(kexp)," ",p.cutoff))
  readkey()
   }



} #{{{main
