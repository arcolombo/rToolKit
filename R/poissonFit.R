#' @title fits a poisson model to determine goodness of fit
#' @description determines if the fit is poisson distributed and fits a poisson regression model
#' @param kexp kallisto experiment at repeat stage level
#' @import PoissonSeq
#' @import pvclust
#' @import dendsort
#' @import ComplexHeatmap
#' @import grid
#' @export 
#' @return data and some images fit to poisson log linear model
poissonFit<-function(kexp, bundleID="tx_id",cutoff=1,comparison="LSC",control="pHSC",p.value=0.05,whichDelta=c("delta1","delta2")){
  whichDelta<-match.arg(whichDelta,c("delta1","delta2"))
##uses poissonSeq to fit the data. 
  kexp<-kexp2Group(kexp,comparison=comparison,control=control)
  kexp<-omitSU583(kexp)
  cnts<-collapseBundles(kexp,"tx_id",read.cutoff=cutoff)
  seq.depth<-PS.Est.Depth(cnts)
  inDat<-list()
  inDat$n<-cnts
  y<-matrix(nrow=ncol(kexp))
  rownames(y)<-colnames(kexp)
  compID<-grep(comparison,colnames(kexp))
  y[compID,]<-2
  contID<-grep(control,colnames(kexp))
  y[contID,]<-1
  y<-as.vector(y)
  inDat$y<-y
  inDat$type<-"twoclass"
  inDat$pair<-FALSE
  inDat$gname<-rownames(cnts)
  res<-PS.Main(inDat) 
###it does class 2 vs class 1 2/1 ########### 

drawDF(kexp,res=res,cutoff=cutoff,fromList=FALSE)
 readkey()
###composition of global identities delta1 LSC-pHSC
  topNames<-as.character(res$gname[which(res$fdr<=0.05)])
  topDE<-res$log.fc[which(res$fdr<=0.05)]
  tt<-data.frame(names=topNames,
               logfc=topDE,
               stringsAsFactors=FALSE)
  rownames(tt)<-topNames
  tt<- tt[!is.infinite(tt$logfc),]
  plotFrequency(kexp,topNames=topNames,topDE=tt,whichDelta=whichDelta)

#### beeswarm and boxplot of global stage
  if(whichDelta=="delta1"){
  plot.new()
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(tt$logfc, main=expression(paste(Delta,"(pHSC,LSC) PoissonSeq")),xlab=paste0("LSC-pHSC p.val",p.value) )
  par(new=TRUE)
  boxplot(tt$logfc,medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
 } else {
  plot.new()
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  beeswarm(tt$logfc, main=expression(paste(Delta,"(LSC,Blast) PoissonSeq")),xlab=paste0("Blast-LSC p.val",p.value) )
  par(new=TRUE)
  boxplot(tt$logfc,medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
   }

####beeswarm by transcript biotypes
biotype_fc<-split(tt,rowRanges(kexp)[rownames(tt)]$tx_biotype)

for(i in 1:length(biotype_fc)){
  yyz<-biotype_fc[[i]]
  familyName<-unique(rowRanges(kexp)[yyz$names]$gene_biotype)
  if(nrow(yyz)>2){
  drawDF(kexp,res=yyz,fromList=TRUE)
  }
  plot.new()
  layout(mat=matrix(c(1,2),ncol=2,byrow=TRUE))
  if(whichDelta=="delta1"){
  beeswarm(yyz$logfc, main=expression(paste(Delta,"(pHSC,LSC) PoissonSeq")),xlab=paste0("LSC-pHSC p.val",p.value," ",familyName) )
  } else {
    beeswarm(yyz$logfc, main=expression(paste(Delta,"(LSC,Blast) PoissonSeq")),xlab=paste0("Blast-LSC p.val",p.value," ",familyName) )
  }
  par(new=TRUE)
  boxplot(yyz$logfc,medcol="red",boxcol="red",whiskcol="red")
  axis(side=4)
  readkey()
  

  }


}
