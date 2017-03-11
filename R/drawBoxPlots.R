#' @title Box plots and custom imaging
#' @description this will plot the box plot of the magnitude of logFC of a DE list and compare medians of the comparison group (one-sided).  this will show the global upregulation of the DE list if the one sided median test is significant.
#' @param kexp this is the non-stage level kexp for 3 stages 
#' @param comparison the comparison group
#' @param control the reference group
#' @param comparison2 if more than one comparison
#' @param numberComparisons  integer, the number of comparisons
#' @param read.cutoff integer, to filter for repeatWise Analysis
#' @param adjustBy either none,BH,BY, or holm for lmFit call
#' @param wilcox.Alternative either less, which is testing the alternative that the comparison is less than the control, or greater which tests the alternative hypothesis that the comparison is greater than control.
#' @param grob.comparison1 this is the text printed to plot that verbose the alternative
#' @param grob.comparison2 this is the second printed explanation of the alternative hypothesis printed on plot testing medians
#' @param testMedians  boolean , whether to include a median test of logFC from the repeat linear modeling.
#' @import arkas
#' @import ggplot2
#' @import grid
#' @import ggthemes
#' @export
drawBoxPlots<-function(kexp,comparison="pHSC",control="LSC",comparison2="Blast",numberComparisons=1,read.cutoff=3,adjustBy="BH",wilcox.Alternative=c("less","greater"),grob.comparison1=paste0("Median(",control,")<Median(",comparison,") p.value: "),grob.comparison2=paste0("Median(",control,")<Median(",comparison2,") p.value: ") ,testMedians=FALSE,title1=NULL,xlab1=NULL,ylab1=NULL,xlab2=NULL  ){
  if(is.null(title1)==TRUE){
   title1<-"Differential Expression Boxplot"
  }
   theme_set(theme_tufte())
  title2<-title1
 
  if(is.null(xlab1)==TRUE){
  xlab1<-"Stage"
  }
  if(is.null(ylab1)==TRUE){
  ylab1<-"|logFC|"
  }
  if(is.null(xlab2)==TRUE){
  xlab2="Stage"
  }
##the density should be absolute value of log FC overlayed across each other. then a wilcoxon sum test should be done. 
  wilcox.Alternative<-match.arg(wilcox.Alternative,c("less","greater"))
  kexp1<-kexp2Group(kexp,comparison=comparison,control=control)
  design<-metadata(kexp1)$design
  gwa1<-repeatWiseAnalysis(kexp1,design=design,adjustBy=adjustBy,read.cutoff=read.cutoff)
  write.csv(gwa1$top,file=paste0("repeatWiseAnalysis.BoxPlot.",comparison,"v.",control,".",adjustBy,".csv"))
  pHSC<-gwa1$top$logFC[which(gwa1$top$logFC>0)]
  lsc<-gwa1$top$logFC[which(gwa1$top$logFC<=0)]
  wilcox.pvalue<-signif(wilcox.test(pHSC,abs(lsc),wilcox.Alternative)$p.value,2)
  lines1=c(rep(comparison,length(pHSC)),rep(control,length(lsc)))
  dat1<-data.frame(logFC=c(abs(pHSC),LSC=abs(lsc)),Stage=lines1)

 dat1$Stage<- factor(dat1$Stage,levels(dat1$Stage)[c(which(levels(dat1$Stage)==comparison),which(levels(dat1$Stage)==control) )])


 if(testMedians==FALSE){
 pp1<-(ggplot(dat1,aes(x=Stage,y=logFC,fill=Stage ))+ggtitle(paste0("Differential Repeat Expression |logFC| ",comparison,"-",control))+scale_fill_manual(values=c("green","red"))+stat_boxplot(aes(Stage,logFC),geom='errorbar',linetype=1,width=0.5)+geom_boxplot(aes(Stage,logFC),outlier.colour=NA))
   pp1<-pp1+geom_point(position=position_jitter(width=0.2),alpha=0.4 )
  print(pp1)
  }else if(testMedians==TRUE){

#  my_grob= grobTree(textGrob(paste0(grob.comparison1,wilcox.pvalue ),x=0.7,y=0.95,gp=gpar(col="black",fontsize=10,fontface="italic")))
  pp1<-(ggplot(dat1,aes(x=Stage,y=logFC,fill=Stage ))+ggtitle(paste0(title1," ",comparison,"-",control), subtitle=paste0(grob.comparison1,wilcox.pvalue   ))+scale_fill_manual(values=c("green","red"))+stat_boxplot(aes(Stage,logFC),geom='errorbar',linetype=1,width=0.5)+geom_boxplot(aes(Stage,logFC),outlier.colour=NA))#+annotation_custom(my_grob,ymax=6,xmax=2) )
pp1<-pp1+xlab(xlab1)+ylab(ylab1)+guides(fill=guide_legend(title="Clonal Stages"))
 pp1<-pp1+geom_point(position=position_jitter(width=0.2),alpha=0.4,show.legend=FALSE )
  print(pp1)
   }
  readkey()
 if(numberComparisons==2){
  kexp2<-kexp2Group(kexp,comparison2,control)
  gwa<-repeatWiseAnalysis(kexp2,design=metadata(kexp2)$design,adjustBy=adjustBy)
  write.csv(gwa$top,file=paste0("repeatWiseAnalysis.Boxplot.",adjustBy,".",comparison2,"v.",control,".csv"))

  blast<-gwa$top$logFC[which(gwa$top$logFC>0)]
  lsc<-gwa$top$logFC[which(gwa$top$logFC<=0)]
  wilcox.pvalue2<-signif(wilcox.test(abs(blast),abs(lsc),wilcox.Alternative)$p.value,2)
  lines=c(rep(comparison2,length(blast)),rep(control,length(lsc)) )
  dat<-data.frame(logFC=c(blast,LSC=abs(lsc)),Stage=lines)
  dat$Stage<- factor(dat$Stage,levels(dat$Stage)[c(which(levels(dat$Stage)==control),which(levels(dat$Stage)==comparison2) )])


  if(testMedians==FALSE){
  pp<-(ggplot(dat,aes(x=Stage,y=logFC,fill=Stage ))+ggtitle(paste0("Differential Repeat Expression |logFC| ",comparison2,"-",control))+scale_fill_manual(values=c("red","lightblue"))+stat_boxplot(aes(Stage,logFC),geom='errorbar',linetype=1,width=0.5)+geom_boxplot(aes(Stage,logFC),outlier.colour=NA)+stat_summary(aes(Stage,logFC),fun.y=mean,geom="point",size=2)+stat_summary(aes(Stage,logFC),fun.data=mean_se,geom="errorbar",width=0.74 ))
   pp<-pp+geom_point(position=position_jitter(width=0.2),alpha=0.4 )
  print(pp)
   }else if(testMedians==TRUE){
   #my_grob= grobTree(textGrob(paste0(grob.comparison2,wilcox.pvalue2 ),x=0.7,y=0.95,gp=gpar(col="black",fontsize=10,fontface="italic")))
  pp<-(ggplot(dat,aes(x=Stage,y=logFC,fill=Stage ))+ggtitle(paste0(title2," ",comparison2,"-",control),subtitle=paste0(grob.comparison2,wilcox.pvalue2 ) )+scale_fill_manual(values=c("red","lightblue"))+stat_boxplot(aes(Stage,logFC),geom='errorbar',linetype=1,width=0.5)+geom_boxplot(aes(Stage,logFC),outlier.colour=NA))# + annotation_custom(my_grob,ymax=6,xmax=2 ) )
  pp<-pp+xlab(xlab2)+ylab(ylab1)+guides(fill=guide_legend(title="Clonal Stages"))

  pp<-pp+geom_point(position=position_jitter(width=0.2),alpha=0.4,show.legend=FALSE ) 
  print(pp)
   }
  
  readkey()

  pdf(paste0("Boxplot_Stage_Plots_",comparison,".",control,".",comparison2,".pdf"))
  print(pp1)
  print(pp)
dev.off()
 }else{
  pdf(paste0("Boxplot_Stage_Plots_",comparison,".",control,".pdf"))
  print(pp1)
  dev.off()
  }

}
