#' @title Density plots and custom imaging
#' @param kexp this is the non-stage level kexp for 3 stages 
#' @import arkas
#' @import ggplot2
#' @export
drawDensityPlots<-function(kexp,comparison="pHSC",control="LSC",comparison2="Blast",numberComparisons=1,read.cutoff=3){


kexp1<-kexp2Group(kexp,comparison=comparison,control=control)
design<-metadata(kexp1)$design
gwa1<-repeatWiseAnalysis(kexp1,design=design,adjustBy="none",read.cutoff=read.cutoff)
 write.csv(gwa1$top,file="repeatWiseAnalysis.pHSC.v.LSC.csv")
 pHSC<-gwa1$top$logFC[which(gwa1$top$logFC>0)]
lsc<-gwa1$top$logFC[which(gwa1$top$logFC<=0)]
lines1=c(rep(comparison,length(pHSC)),rep(control,length(lsc)))
dat1<-data.frame(logFC=c(pHSC,LSC=lsc),Stage=lines1)

 print(ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison,"-",control))+scale_fill_manual(values=c("red","green")))
  readkey()




 if(numberComparisons==2){


 kexp2<-kexp2Group(kexp,comparison2,control)
 gwa<-repeatWiseAnalysis(kexp2,design=metadata(kexp2)$design,adjustBy="none")
 write.csv(gwa$top,file="repeatWiseAnalysis.Blast.v.LSC.csv")

 blast<-gwa$top$logFC[which(gwa$top$logFC>0)]
 lsc<-gwa$top$logFC[which(gwa$top$logFC<=0)]
 lines=c(rep(comparison2,length(blast)),rep(control,length(lsc)))
 dat<-data.frame(logFC=c(blast,LSC=lsc),Stage=lines)

  print(ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison2,"-",control))+scale_fill_manual(values=c("blue","red")))
  readkey()

  pdf("Density_Stage_Plots.pdf")
  print(ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison,"-",control))+scale_fill_manual(values=c("red","green")))
  print(ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison2,"-",control))+scale_fill_manual(values=c("blue","red")))
dev.off()
 }else{
  pdf("Density_Stage_Plots.pdf")
  print(ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison,"-",control))+scale_fill_manual(values=c("red","green")))
 }

}
