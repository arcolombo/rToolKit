#' @title this creates an expression time plot of match patients
#' @description for non-paired group factor ANOVA testing, this plots the expression across two groups units of TPM
#' @param kexp this is a kexp for 2 gropu non-paired patient wise, with 2 groups
#' @import ComplexHeatmap
#' @import ggplot2
#' @import ggthemes
#' @export
#' @return images and a pdf/jpeg
geneTest<-function(kexp,geneName="EVI",how=c("cpm","tpm"),saveOut=FALSE,numberPairs=c(1,2,3),read.cutoff=1){

 theme_set(theme_tufte())
 rexp<-kexp
 #rexp<-kexpByTrio(rexp) ##take all triples
 pairs<-colnames(rexp)
 pairs<-unique(sapply(strsplit(pairs,"_"),function(x) x[2]))
 samples<-unique(sapply(strsplit(colnames(rexp),"_"),function(x) x[1]))

 n.samples<-factor(sapply(strsplit(colnames(rexp),"_"),function(x) x[1]))
 number.samples<-sapply(samples,function(x) length(grep(x[1],n.samples)))
 samples<-paste0(samples,"_Risk")

 if(how=="tpm"){
 tpm<-collapseTpm(rexp,"gene_name",minTPM=read.cutoff) ##36 repeat classes collapsed
 ##now for each repeat class and for each pairs, plot a connected dot plot
 tpm<-tpm[which(rownames(tpm)==geneName),]
 stopifnot(nrow(tpm)>0)
 }else{
  dge<-list()
 bundledCounts <- collapseBundles(kexp, bundleID="gene_name",
                                   read.cutoff=1)
  dge <- DGEList(counts=bundledCounts)
  dge <- calcNormFactors(dge)
  tpm<-cpm(dge)
  tpm<-tpm[which(rownames(tpm)==geneName),]
 }

 if(length(tpm)==0){
  return(0)
  }

 # write.csv(tpm,file=paste0(geneName,"_Gene-Matched_Triplicates.TPM.PatientTrioPlots.csv"))
  pair.tpm<-t(tpm)
  pair.tpm<-log2(1+pair.tpm)  
  

 m<-as.matrix(pair.tpm)
 group<-data.frame() 
for(i in 1:length(samples)){
 group<-c(group,rep(samples[i],number.samples[i]))
 }
 
  group<-unlist(group)
  data=data.frame(y=t(m),group=factor(group,levels=samples))

###add repeated measures
#for(i in 1:length(pairs)){


#}
#######

  fit<-lm(y~group,data)
 ANOVA<-anova(fit)
  anova.pvalue<-ANOVA[1,5]
 pair.T<-pairwise.t.test(data$y,data$group,p.adj='bonferroni')
 if(numberPairs==3){
 pair.pHSC.LSC<-pair.T$p.value[1,1]
 pair.Blast.LSC<-pair.T$p.value[2,2]
 pair.Blast.pHSC<-pair.T$p.value[2,1]
   pairWise.DF<-data.frame(ANOVA=anova.pvalue,pHSC.v.LSC=pair.pHSC.LSC,Blast.v.LSC=pair.Blast.LSC,Blast.v.pHSC=pair.Blast.pHSC) 
 }else if(numberPairs==1){
   pairWise.DF<-data.frame(ANOVA=anova.pvalue,Pair=pair.T$p.value)
 }
 names(pairWise.DF)<-gsub(".v.","-",names(pairWise.DF))
  alpha<-which(pairWise.DF<=0.05)
  pchLabels<-as.numeric(rep(0,ncol(pairWise.DF)))
  pchLabels[alpha]<-8
 #names(pairWise.DF)[alpha]<-paste0(names(pairWise.DF)[alpha],"**")
print(pairWise.DF)
 if(saveOut==TRUE){
  write.csv(pairWise.DF,file=paste0(geneName,"_ANOVA.csv"))
 }
 readkey() 
 
  pp1<-ggplot(data,aes(x=group,y=y,file=group))+ggtitle(paste0(geneName," Coding Gene Expression"))+stat_boxplot(aes(group,y),geom='errorbar',linetype=1,width=0.5)+geom_boxplot(aes(group,y))+theme(text=element_text(size=20))+xlab("Low/High Risk MDS")+ylab(paste0(geneName," Expression Level (Log2 ",toupper(how),")"))
  pp1<-pp1+geom_point(position=position_jitter(width=0.2),alpha=0.4)
 print(pp1)
 readkey()
   gn.mt<-t(m)
  colnames(gn.mt)<-geneName
 geneHeat<-(Heatmap(gn.mt,name=paste0("Log(1+",toupper(how),")")))
 print(geneHeat)
 readkey()
  if(saveOut==TRUE){
 
   pdf(paste0(geneName,"_ANOVA.pdf"))
   print(pp1)
   print(geneHeat)
   dev.off()
 }
 cat("done.")
}##main
