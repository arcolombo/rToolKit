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
  pair.tpm<-t(tpm)
  
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
     par(xaxt="n")
    matplot(DF, type = c("b"),pch=1,col = 1:ncol(DF),main=paste0(cols," TPM Repeat Expression"),ylab="TPM" ) #plot
   legend("topright",legend=colnames(DF),col=1:ncol(DF),pch=1)
  axis(1,at=seq(1,3,1),labels=FALSE)
  lablist.x<-c(stage1,stage2,stage3)
 text(x=seq(1,3,1),par("usr")[3]-0.2,labels=lablist.x,xpd=TRUE,srt=45,pos=1)
  readkey()
   leadTitle<-gsub("/","",cols)
   leadTitle<-gsub(" ","_",leadTitle)
   if(printWhat=="pdf"){
   pdf(paste0(leadTitle,"patientTrio_RepeatPlot.pdf"))
   par(xaxt="n")
   matplot(DF, type = c("b"),pch=1,col = 1:ncol(DF),main=paste0(cols," TPM Repeat Expression"),ylab="TPM" ) #plot
   legend("topright",legend=colnames(DF),col=1:ncol(DF),pch=1)
   axis(1,at=seq(1,3,1),labels=FALSE)
   lablist.x<-c(stage1,stage2,stage3)
   text(x=seq(1,3,1),par("usr")[3]-0.2,labels=lablist.x,xpd=TRUE,srt=45,pos=1)
  dev.off()
  }else{
   jpeg(paste0(leadTitle,"patientTrio_RepeatPlot.jpeg"))
   par(xaxt="n")
   matplot(DF, type = c("b"),pch=1,col = 1:ncol(DF),main=paste0(cols," TPM Repeat Expression"),ylab="TPM" ) #plot
   legend("topright",legend=colnames(DF),col=1:ncol(DF),pch=1)
   axis(1,at=seq(1,3,1),labels=FALSE)
   lablist.x<-c(stage1,stage2,stage3)
   text(x=seq(1,3,1),par("usr")[3]-0.2,labels=lablist.x,xpd=TRUE,srt=45,pos=1)
  dev.off()
    }##Jpeg
  } ##tx biotypes
}##main
