#' @title Density plots and custom imaging
#' @param kexp this is the non-stage level kexp for 3 stages 
#' @import arkas
#' @import ggplot2
#' @import grid
#' @export
drawDensityPlots<-function(kexp,comparison="pHSC",control="LSC",comparison2="Blast",numberComparisons=1,read.cutoff=3,adjustBy="BH",wilcox.Alternative=c("less","greater"),grob.comparison1=c("Md(High Grade)<Md(Low Grade) p.value: "),grob.comparison2=c("Md(High Grade)<Md(Low Grade p.value: ") ,testMedians=FALSE  ){

##the density should be absolute value of log FC overlayed across each other. then a wilcoxon sum test should be done. 
  wilcox.Alternative<-match.arg(wilcox.Alternative,c("less","greater"))
kexp1<-kexp2Group(kexp,comparison=comparison,control=control)
design<-metadata(kexp1)$design
gwa1<-repeatWiseAnalysis(kexp1,design=design,adjustBy=adjustBy,read.cutoff=read.cutoff)
 write.csv(gwa1$top,file=paste0("repeatWiseAnalysis.",comparison,"v.",control,".",adjustBy,".csv"))
 pHSC<-gwa1$top$logFC[which(gwa1$top$logFC>0)]
lsc<-gwa1$top$logFC[which(gwa1$top$logFC<=0)]
wilcox.pvalue<-signif(wilcox.test(abs(pHSC),abs(lsc),wilcox.Alternative)$p.value,2)
lines1=c(rep(comparison,length(pHSC)),rep(control,length(lsc)))
dat1<-data.frame(logFC=c(abs(pHSC),LSC=abs(lsc)),Stage=lines1)

 if(testMedians==FALSE){
 pp1<-(ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison,"-",control," FDR: ",adjustBy))+scale_fill_manual(values=c("red","green")))
  print(pp1)
  }else if(testMedians==TRUE){

  my_grob= grobTree(textGrob(paste0(grob.comparison1,wilcox.pvalue ),x=0.7,y=0.95,gp=gpar(col="black",fontsize=10,fontface="italic")))
 pp1<-(ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison,"-",control," FDR: ",adjustBy))+scale_fill_manual(values=c("red","green"))+annotation_custom(my_grob) )
 print(pp1)
   }
  readkey()
 if(numberComparisons==2){
  kexp2<-kexp2Group(kexp,comparison2,control)
  gwa<-repeatWiseAnalysis(kexp2,design=metadata(kexp2)$design,adjustBy=adjustBy)
  write.csv(gwa$top,file=paste0("repeatWiseAnalysis.",adjustBy,".",comparison2,"v.",control,".csv"))

  blast<-gwa$top$logFC[which(gwa$top$logFC>0)]
  lsc<-gwa$top$logFC[which(gwa$top$logFC<=0)]
  wilcox.pvalue2<-signif(wilcox.test(abs(blast),abs(lsc),wilcox.Alternative)$p.value,2)
  lines=c(rep(comparison2,length(blast)),rep(control,length(lsc)))
  dat<-data.frame(logFC=c(blast,LSC=abs(lsc)),Stage=lines)
  if(testMedians==FALSE){
  pp<-(ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison2,"-",control," FDR: ",adjustBy))+scale_fill_manual(values=c("blue","red")))
  print(pp)
   }else if(testMedians==TRUE){
   my_grob= grobTree(textGrob(paste0(grob.comparison2,wilcox.pvalue ),x=0.7,y=0.95,gp=gpar(col="black",fontsize=10,fontface="italic")))
 pp<-(ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle(paste0("Differential Repeat Expression ",comparison2,"-",control," FDR: ",adjustBy))+scale_fill_manual(values=c("blue","red"))+annotation_custom(my_grob) )
   print(pp)
   }
  
  readkey()

  pdf("Density_Stage_Plots.pdf")
  print(pp1)
  print(pp)
dev.off()
 }else{
  pdf("Density_Stage_Plots.pdf")
  print(pp1)
  dev.off()
  }

}
