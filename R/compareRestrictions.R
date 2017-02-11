#' @title compareRestrictions analysis of qusage
#' @description compares the qusage results from two different experiments and finds the correlation of means for each pathway level after standardizing the means.
#' @param kexp1 this is the first experiment you wish to analyze
#' @param kexp2 the second experiment you wish to examine across.  
#' @param geneSetPath the path to the gmt files which contain the geneSet.gmt files
#' @param compare1 the first comparison group in kexp1
#' @param compare2 the second comparison group in kexp2
#' @param controls1 the first control/reference group in kexp1
#' @param controls2 the second control/reference group in kexp2
#' @param paired1 boolean indicating if kexp1 is paired samples
#' @param paired2 boolean indicating if kexp2 is paired samples
#' @param read.cutoff min threshold of count data
#' @param MsigDB character of the .gmt file to load
#' @param saveToFile boolean will save a pdf and write the data to .csv
#' @import qusage
#' @import ggplot2
#' @import ggpmisc
#' @import ggthemes
#' @import arkas
#' @export
compareRestrictions<-function(kexp1,kexp2,geneSetPath="~/Documents/Arkas-Paper-Data/MSigDB/MsigDb_all/",compare1="Low",compare2="pHSC",controls1="High",controls2="LSC",paired1=FALSE,paired2=TRUE,read.cutoff=2,MsigDB=NULL,saveToFile=FALSE){

 geneSet<-MsigDB
 comparison1<-compare1
 comparison2<-compare2
   kexp1<-kexp2Group(kexp1,comparison=compare1,control=controls1)
    kexp2<-kexp2Group(kexp2,comparison=compare2,control=controls2)
  qusageMean1<-qusageMean(kexp=kexp1,read.cutoff=read.cutoff,comparison=comparison1,controls=controls1,geneSet=geneSet,paired=paired1,geneSetPath=geneSetPath)
   qusageMean2<-qusageMean(kexp=kexp2,read.cutoff=read.cutoff,comparison=comparison2,controls=controls2,geneSet=geneSet,paired=paired2,geneSetPath=geneSetPath)

 ##standardize for correlation tests is not needed!
 #  norm1<- lapply(qusageMean1,function(x) (x-mean(x))/sd(x))
 # norm2<-lapply(qusageMean2,function(x) (x-mean(x))/sd(x))
 norm1<-qusageMean1
 norm2<-qusageMean2

 corVec<-matrix(nrow=length(names(norm1)),ncol=2)
 rownames(corVec)<-names(norm1)
 colnames(corVec)<-c("Correlation","p.value")
# if(saveToFile==TRUE){
# write.csv(corVec,file="MDS-AML-RestrictionFactor-Analysis-Correlations-pvalues.csv")
# }
 ##match intersect
  for(i in 1:length(names(norm1))){
  n1.id<-intersect(names(norm1[[i]]),names(norm2[[i]]))
  n1<-norm1[[i]][n1.id]
  n2<-norm2[[i]][n1.id]
  corVec[i,1]<-cor(n1,n2)
  corVec[i,2]<-cor.test(n1,n2)$p.value
 df<-data.frame(MDS=(n1),AML=(n2))
 #my.formula <- y ~ x
  p<-  ggplot(df, aes(x = MDS, y = AML)) +
  geom_point(shape = 19, size = 2 ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm',formula='y~x',se=FALSE) +
  stat_poly_eq(aes(label=paste(..rr.label..)),
                   label.x.npc="right",label.y.npc=1.15,
                   formula='y~x',parse=TRUE,size=5)+
  stat_fit_glance(method='lm',
                  method.args=list(formula='y~x'),
                  geom='text',
                  aes(label=paste("P.value = ",signif(..p.value..,digits=4),sep="")),
   label.x.npc='right',
   label.y.npc=1.35,size=5)+
  ggtitle(paste0("MDS-AML Restriction Factor Correlations: ",names(norm1)[i])) +
  xlab("MDS") +
  ylab("AML") +
  scale_colour_tableau("tableau10") +
  #geom_text(x = -1, y = 0.5,
   #         label = corr_eqn(df$MDS,
    #                         df$AML), parse = TRUE) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(lineheight = .8, face = "bold", vjust = 1),
        axis.text.x = element_text(size = 11, vjust = 0.5,
                                   hjust = 1, colour = 'black'),
        axis.text.y = element_text(size = 11, colour = 'black'),
        axis.title = element_text(size = 10, face = 'bold'),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(colour = 'black', size = 1),
        panel.background = element_blank())

 print(p)
 readkey()
 if(saveToFile==TRUE){
pdf(paste0("MDS-AML_ImmuneComplexes_pHSC-Blast-Restriction_Factor_Correlations_",gsub(" ","",names(norm1)[i]),".pdf" ))
# print(p)
 print(p)
 dev.off() 
 }
}


 if(saveToFile==TRUE){
 write.csv(corVec,file="MDS-AML-RestrictionFactor-Analysis-Correlations-pvalues.csv")
 }


} #main 

#' @title qusage mean activity
#' @description calculates the mean pathway activity for a gene set
#' @import edgeR
#' @import qusage
qusageMean<-function(kexp=NULL,read.cutoff=read.cutoff,comparison=NULL,controls=NULL,geneSet=geneSet,paired=FALSE,geneSetPath=geneSetPath  ){


  counts<-collapseBundles(kexp,"gene_name",read.cutoff=read.cutoff)
  ##BSN
  counts<-counts[,which(colSums(counts)>1)]
   dge<-DGEList(counts=counts)
   dge<-calcNormFactors(dge)
   expr<-cpm(dge,normalized.lib.sizes=TRUE,log=FALSE)
  expr2<-log2(1+expr)

  cN<-colnames(expr2)
  cN<-strsplit(cN,"_")
  labels<-unlist(lapply(cN,function(x) x[1]))
 
   contrast<-paste0(comparison,"-",controls)
   geneSets<-read.gmt(paste0(geneSetPath,geneSet))
  ##alphabetical order the geneSets must be
    geneSets<-geneSets[sort(names(geneSets))]
  
 if(paired==TRUE){
   pairs<-unlist(lapply(cN,function(x) x[2]))
   pairs.id<-match(toupper(pairs),toupper(pairs))
   qs.results<-qusage(expr2,labels,contrast,geneSets,pairVector=pairs.id)
  }else if(paired==FALSE){
   qs.results<-qusage(expr2,labels,contrast,geneSets)
  }
  tt<-qs.results$pathways
 out<-lapply(tt,function(x) sort((qs.results$mean[x])) )
 return(out)
}

