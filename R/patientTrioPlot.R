#' @title this creates an expression time plot of match patients
#' @description for paired time series, this plots the expression across time points in units of TPM. This plots each repeat biotype and every patient trio across stages.
#' @param kexp this is a kexp for three stages patient wise, with 3 time-stages
#' @param stage1 character pHSC
#' @param stage2 character LSC
#' @param stage3 character Blast
#' @export
#' @return images and a pdf/jpeg
patientTrioPlot<-function(kexp,stage1="pHSC",stage2="LSC",stage3="Blast",printWhat="jpeg"){

 printWhat=match.arg(printWhat,c("pdf","jpeg"))
 rexp<-findRepeats(kexp)
 rexp<-kexpByTrio(rexp) ##take all triples
 pairs<-colnames(rexp)
 pairs<-unique(sapply(strsplit(pairs,"_"),function(x) x[2]))
 tpm<-collapseTpm(rexp,"tx_biotype") ##36 repeat classes collapsed
 ##now for each repeat class and for each pairs, plot a connected dot plot
  write.csv(tpm,file="TxBiotype-Repeat-Matched_Triplicates.TPM.PatientTrioPlots.csv")
  pair.tpm<-t(tpm)
  pair.tpm<-log2(1+pair.tpm)  
  for(cols in colnames(pair.tpm)){ 
  #grab first
  paird<-pair.tpm[grep(pairs[1],rownames(pair.tpm)),grep(cols,colnames(pair.tpm))]
  DF<-data.frame(paird)
  colnames(DF)<-pairs[1]
  rownames(DF)<-c(stage1,stage2,stage3)
  for(j in 2:length(pairs)){
   paird<-pair.tpm[grep(pairs[j],rownames(pair.tpm)),grep(cols,colnames(pair.tpm))]
  df<-data.frame(paird)
  colnames(df)<-pairs[j]
  rownames(df)<-c(stage1,stage2,stage3)
  DF<-cbind(DF,df)
  }
 m<-as.matrix(DF)
y<-c(as.vector(m[1,]),as.vector(m[2,]),as.vector(m[3,]))
  group<-c(rep(1,7),rep(2,7),rep(3,7))
  data=data.frame(y=y,group=factor(group))
  fit<-lm(y~group,data)
 ANOVA<-anova(fit)
  anova.pvalue<-ANOVA[1,5]
 pair.T<-pairwise.t.test(data$y,data$group,p.adj='bonferroni')
 pair.pHSC.LSC<-pair.T$p.value[1,1]
 pair.Blast.LSC<-pair.T$p.value[2,2]
 pair.Blast.pHSC<-pair.T$p.value[2,1]
  pairWise.DF<-data.frame(ANOVA=anova.pvalue,pHSC.v.LSC=pair.pHSC.LSC,Blast.v.LSC=pair.Blast.LSC,Blast.v.pHSC=pair.Blast.pHSC) 
  x11(width=8,height=8)
     par(mar=par()$mar+c(0,0,0,5),family='Helvetica',xpd=TRUE,cex.main=1)
    matplot(DF, type = c("b"),pch=1,col = c("black","green","orange","red","purple","blue","magenta")  ,xaxt='n',xlab=NA,main=paste0("'",cols,"' Transposable Element Family ANOVA Across Clonal Stages"),ylab="Transposable Element Expression Level (Log TPM)" ) #plot
   legend("topright",legend=colnames(DF),col=c("black","green","orange","red","purple","blue","magenta") ,pch=0.8,title="Patient ID",inset=c(-0.2,0))
  axis(1,at=seq(1,3,1),labels=c(stage1,stage2,stage3))
  lablist.x<-c(stage1,stage2,stage3)
 # text(x=seq(1,3,1),(par("mgp")[1]-3.4),labels=lablist.x,xpd=TRUE,srt=45,pos=1,cex=1.2)

 ### anova stats
 legend("top",legend=paste0(colnames(pairWise.DF),": ",signif(pairWise.DF,1)),pch=1,title="Adj.P.Values" )
 
 



 ###

  readkey()
   leadTitle<-gsub("/","",cols)
   leadTitle<-gsub(" ","_",leadTitle)
   if(printWhat=="pdf"){
   pdf(paste0(leadTitle,"patientTrio_RepeatPlot.pdf"),family='Helvetica',height=8,width=9)
    
  par(mar=par()$mar+c(0,0,0,5),family='Helvetica',xpd=TRUE,cex.main=1)
    matplot(DF, type = c("b"),pch=1,col = c("black","green","orange","red","purple","blue","magenta")  ,xaxt='n',xlab=NA,main=paste0("'",cols,"' Transposable Element Family ANOVA Across Clonal Stages"),ylab="Transposable Element Expression Level (Log TPM)" ) #plot
   legend("topright",legend=colnames(DF),col=c("black","green","orange","red","purple","blue","magenta") ,pch=0.8,title="Patient ID",inset=c(-0.2,0))
  axis(1,at=seq(1,3,1),labels=c(stage1,stage2,stage3))
 
  ### anova stats
 legend("top",legend=paste0(colnames(pairWise.DF),": ",signif(pairWise.DF,1)),pch=1,title="Adj.P.Values" )
  dev.off()
  }else{
   jpeg(paste0(leadTitle,"patientTrio_RepeatPlot.jpeg"))
 
 par(mar=par()$mar+c(0,0,0,5),family='Helvetica',xpd=TRUE,cex.main=1)
    matplot(DF, type = c("b"),pch=1,col = c("black","green","orange","red","purple","blue","magenta")  ,xaxt='n',xlab=NA,main=paste0("'",cols,"' Transposable Element Family ANOVA Across Clonal Stages"),ylab="Transposable Element Expression Level (Log TPM)" ) #plot
   legend("topright",legend=colnames(DF),col=c("black","green","orange","red","purple","blue","magenta") ,pch=0.8,title="Patient ID",inset=c(-0.2,0))
  axis(1,at=seq(1,3,1),labels=c(stage1,stage2,stage3))
  lablist.x<-c(stage1,stage2,stage3)
 # text(x=seq(1,3,1),(par("mgp")[1]-3.4),labels=lablist.x,xpd=TRUE,srt=45,pos=1,cex=1.2)

 ### anova stats
 legend("top",legend=paste0(colnames(pairWise.DF),": ",signif(pairWise.DF,1)),pch=1,title="Adj.P.Values" )


  dev.off()
    }##Jpeg
  } ##tx biotypes
  cat("done.\n")
}##main
