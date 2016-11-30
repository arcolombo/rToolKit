#' @title Density plots and custom imaging
#' @param kexp this is the non-stage level kexp for 3 stages 
drawDensityPlots<-function(kexp){


kexp1<-kexp2Group(kexp,"pHSC","LSC")
gwa1<-repeatWiseAnalysis(kexp,design=metadata(kexp)$design,adjustBy="none")
 pHSC<-gwa$top$logFC[which(gwa$top$logFC>0)]
lsc<-gwa$top$logFC[which(gwa$top$logFC<=0)]
lines1=c(rep("pHSC",length(blast)),rep("LSC",length(lsc)))
dat1<-data.frame(logFC=c(pHSC,LSC=lsc),Stage=lines1)

ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle("Differential Repeat Expression pHSC-LSC")+scale_fill_manual(values=c("red","green"))
 readkey()







kexp2<-kexp2Group(kexp,"Blast","LSC")
gwa<-repeatWiseAnalysis(kexp,design=metadata(kexp)$design,adjustBy="none")
 blast<-gwa$top$logFC[which(gwa$top$logFC>0)]
lsc<-gwa$top$logFC[which(gwa$top$logFC<=0)]
lines=c(rep("Blast",length(blast)),rep("LSC",length(lsc)))
dat<-data.frame(logFC=c(blast,LSC=lsc),Stage=lines)

ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle("Differential Repeat Expression Blast-LSC")+scale_fill_manual(values=c("blue","red"))
 readkey()

pdf("Density_Stage_Plots.pdf")
ggplot(dat1,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle("Differential Repeat Expression pHSC-LSC")+scale_fill_manual(values=c("red","green"))
ggplot(dat,aes(x=logFC,fill=Stage ))+geom_density(alpha=0.5)+ggtitle("Differential Repeat Expression Blast-LSC")+scale_fill_manual(values=c("blue","red"))
dev.off()


}
