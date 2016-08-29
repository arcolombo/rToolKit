#' @title find the repeats falling into orders
#' @description  find the repeats fallingo into clonal ordering with the alternative model and hypothesis model input for element selection.  
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
findModel<-function(kexp,globalMax="pHSC", globalMin="LSC",patientID=NULL,repeatType=NULL,read.cutoff=1,outputReport=FALSE,outputDir=NULL){
 
 #the input is a clonal kexp, findModel will split
  if(is.null(patientID)==FALSE){
 kexp<-kexpByPatient(kexp,patientID=patientID)
  }
  globalmax<-grep(globalMax,colnames(kexp))
  globalmin<-grep(globalMin,colnames(kexp))
   midID<-which(colnames(kexp)!=colnames(kexp)[globalmax])
  globalmid<-midID[which(colnames(kexp)[midID] !=colnames(kexp)[globalmin])]
   globalMid<-colnames(kexp)[globalmid]
  
  mt<-collapseBundles(kexp,"tx_id",read.cutoff=read.cutoff)
  ##quety repeats with globalMax > globalmid> globalMin
  id1<- which(mt[,globalmax]>mt[,globalmin])
   tt<-mt[id1,]
   stopifnot(all(tt[,globalmax]>tt[,globalmin])==TRUE)
  id2<-which(tt[,globalmax]>tt[,globalmid])
  t2<-tt[id2,]
   stopifnot(all(t2[,globalmax]>t2[,globalmid]))
 
  ##fix me add all query groups
  ## globalMax =globalmid>globalMin
  ## gobalMax>mid=min
  ##Mid > max > min
  ## max = mid = min
 ## Max > min > mid
 ## Max < min < mid

 ##inverted
  ##min > max > mid
  ##min = mid > max
  ## 

  hypothesis<-t2
  ##find annotations
  df<-data.frame(names=rownames(hypothesis),
                 tx_biotype=rowRanges(kexp)[rownames(hypothesis)]$tx_biotype,
                 gene_biotype=rowRanges(kexp)[rownames(hypothesis)]$gene_biotype)
 
  df2<-cbind(hypothesis,df)

  dTab<-(table(df[,2]))
  #####freq bar chart of nominal category variables
  barplot(dTab,main=paste0("TxBiotypes N=",nrow(hypothesis),"/",nrow(mt)," ",colnames(kexp)[globalmax],">",globalMid,">",colnames(kexp)[globalmin],">",read.cutoff ),las=2,cex.names=0.6)
  readkey()

  ######## relative freq bar chart of category variables
   barplot((dTab/sum(dTab))*100, main=paste0("TxBio RF N=",nrow(hypothesis),"/",nrow(mt)," ",colnames(kexp)[globalmax],">",globalMid,">",colnames(kexp)[globalmin],">",read.cutoff ),las=2,cex.names=0.6)
  readkey()


  ######gene biotype
  barplot(table(df[,3]),main=paste0("Genetypes Freq N=",nrow(hypothesis),"/",nrow(mt)," ",colnames(kexp)[globalmax],">=",globalMid,">",colnames(kexp)[globalmin],">0"),las=2,cex.names=0.6)


  #########pie chart 
  pie(dTab[which(dTab>read.cutoff)] )
  dTabFreq <- prop.table(dTab[which(dTab>1)] )
  textRad  <- 0.5
  angles   <- dTabFreq * 2 * pi
  csAngles <- cumsum(angles)
  csAngles <- csAngles - angles/2
  textX    <- textRad * cos(csAngles)
  textY    <- textRad * sin(csAngles)
  text(x=textX, y=textY, labels=round(dTabFreq,3))
  title(paste0("TxBiotypes N=",nrow(hypothesis),"/",nrow(mt)," ",colnames(kexp)[globalmax],">=",globalMid,">",colnames(kexp)[globalmin],">",read.cutoff ))
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
  title(paste0("GeneBiotypes N=",nrow(hypothesis),"/",nrow(mt)," ",colnames(kexp)[globalmax],">=",globalMid,">",colnames(kexp)[globalmin],">0"))
  readkey()


  boxplot(asinh(df2[,c(globalmax,globalmid,globalmin)]),medcol="red",boxcol="red",whiskcol="red")
   axis(side=4)
  par(new=TRUE)
  beeswarm(asinh(df2[,c(globalmax,globalmid,globalmin)]), pch=16,xlab=paste0(colnames(kexp)[globalmax],">",globalMid,">",colnames(kexp)[globalmin],">1") ) 
    readkey()
  
   
   ##FIX ME:  print out the top 40-50 repeats
 ##print out

 if(outputReport==TRUE){
 write.csv(df2,file=paste0(outDir,"/",globalMax,"_",globalMid,"_",globalMin,".csv",row.names=TRUE))
  ##FIX ME:  create a pdf with every repeat 1 pdf

 } else {
  return(df2) 
 }



  ##FIX ME:  should not the groups  A > B > C   and A < B < C  be equal to N? 
 ###FIX ME:  do a multiplot of pHSC > LSC and LSC > pHSC  show all combos
  ###FIX ME ::  N(pHSC>LSC) + N(LSC>pHSC ) = N_total
} #{{{main
