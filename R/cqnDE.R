#' @title CQN and edgeR analysis
#' @description this runs CQN and edgeR on a kallistoExperiment with dynamic contrasts to select  CQN normalization across clonal stages.  This must have replicates so we are comparing samples and not specific patients.
#' @param kexp a kallistoExperiment preferably only repeats, also the repeat kexp must not be TMM Normalized
#' @param contrastMatrix a matrix with contrasts, this is specific and this function is created without fail safes
#' @param design a matrix to compare, 2 group treatment factor
#' @param comparison  the column group either pHSC , LSC ,or Blast
#' @param control  not comparison
#' @param group3   not comparison and not control
#' @param cutoff   integer  a minimum cutoff threshold
#' @import cqn
#' @import edgeR 
#' @import quantreg 
#' @importFrom scales alpha
#' @importFrom graphics par grid
#' @importFrom grid unit gpar
#' @import ComplexHeatmap
#' @export
#' @return returns 
cqnDE<-function(kexp,contrastMatrix=NULL,design=NULL,patientID=NULL, comparison=comparison,control=control,group3=NULL,cutoff=2){

  if(is.null(design)==TRUE || is.null(metadata(kexp)$design)==TRUE ) {
  kexp<-kexp2Group(kexp,comparison=comparison,control=control)
   design<-metadata(kexp)$design
   } 
    
 
  message("note the kexp must not be TMM normalized")
  cnts<-collapseBundles(kexp,bundleID="tx_id",read.cutoff=cutoff) 
  dge<-DGEList(counts=cnts)
  sizeFactors.subset<-dge$samples$lib.size #grabs sequence depth
  names(sizeFactors.subset)<-rownames(dge$samples)
  uCovar<-data.frame(length=rowRanges(kexp)[rownames(cnts)]$tx_length,gccontent=rowRanges(kexp)[rownames(cnts)]$gc_content)
  rownames(uCovar)<-rownames(cnts)
  

  stopifnot(all(rownames(cnts)==rownames(uCovar)))
  stopifnot(all(colnames(cnts)==names(sizeFactors.subset)))
  ##FIX ME: use CPM values and not RPKM??
  cqn.subset<-cqn(round(cnts),lengths=uCovar$length,x=uCovar$gccontent,sizeFactors=sizeFactors.subset,verbose=TRUE,eff_len=eff_len,byWhat="tpm")


  par(mfrow=c(1,2))
 cqnplot(cqn.subset, n = 1,lty=1,xlab=(paste0("GC content Samples RPKM ",comparison," ",control)))
 cqnplot(cqn.subset, n = 2, xlab =paste0("length Samples RPKM ",comparison," ",control), lty = 1)
  readkey()
  RPKM.cqn<-cqn.subset$y + cqn.subset$offset #on log2 scale

  RPM <- sweep(log2(cnts + 1), 2, log2(sizeFactors.subset/10^6)) #CPM
  RPKM.std <- sweep(RPM, 1, log2(uCovar$length / 10^3)) #standard RPKM
  
  
  #need to identify groups
  grp1<-pData(kexp)[grepl(comparison,pData(kexp)$ID),]
  grp2<-pData(kexp)[!grepl(comparison,pData(kexp)$ID), ]
  
  whGenes <- which(rowMeans(RPKM.std) >= 2 & uCovar$length >= 100)
  M.std <- rowMeans(as.data.frame(RPKM.std[whGenes, grp1])) - rowMeans(as.data.frame(RPKM.std[whGenes, grp2])) 
  A.std <- rowMeans(RPKM.std[whGenes,])
  M.cqn <- rowMeans(as.data.frame(RPKM.cqn[whGenes, grp1])) - rowMeans(as.data.frame(RPKM.cqn[whGenes, grp2])) 
  A.cqn <- rowMeans(RPKM.cqn[whGenes,])

par(mfrow = c(1,2))
 plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
 main = "Standard RPKM", ylim = c(-4,4), xlim = c(0,12),
 col = alpha("black", 0.55))
 plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
 main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12),
 col = alpha("black", 0.55))
  readkey()


par(mfrow = c(1,2))
 gccontent <- uCovar$gccontent[whGenes]
 whHigh <- which(gccontent > quantile(gccontent, 0.9))
 whLow <- which(gccontent < quantile(gccontent, 0.1))
 plot(A.std[whHigh], M.std[whHigh], cex = 0.55, pch = 16, xlab = "A",
 ylab = "M", main = "Standard RPKM",
 ylim = c(-4,4), xlim = c(0,13), col = "red")
 points(A.std[whLow], M.std[whLow], cex = 0.55, pch = 16, col = "blue")
 plot(A.cqn[whHigh], M.cqn[whHigh], cex = 0.55, pch = 16, xlab = "A",
 ylab = "M", main = "CQN normalized RPKM",
 ylim = c(-6,6), xlim = c(0,19), col = "red")
 points(A.cqn[whLow], M.cqn[whLow], cex = 0.55, pch = 16, col = "blue")
 readkey()

########edgeR differential expression analysis
 
d.mont <- DGEList(counts = cnts, lib.size = sizeFactors.subset, group = sapply(strsplit(pData(kexp)$ID,"_"),function(x) x[1])
, genes = uCovar)

  ##FIX ME add multi group
  d.mont$offset<-cqn.subset$glm.offset
  d.mont.cqn<-estimateGLMCommonDisp(d.mont,design=design)

  efit.cqn<-glmFit(d.mont.cqn,design=design)
  elrt.cqn<-glmLRT(efit.cqn,coef=2)
  tagd<-topTags(elrt.cqn,n=100,adjust="BH",p=0.05)

  #####heatmap ###########
  drawHeatmap(kexp,tags=tagd,byType="counts",cutoff=2)
 

  ###heatmap of repeats
  cnts<-cnts[!grepl("^ERCC",rownames(cnts)),]
  cnts<-cnts[!grepl("^ENST",rownames(cnts)),]
  d.mont <- DGEList(counts = cnts, lib.size = sizeFactors.subset, group = sapply(strsplit(pData(kexp)$ID,"_"),function(x) x[1])
, genes = uCovar[rownames(uCovar)%in%rownames(cnts),])

  ##FIX ME add multi group
  d.mont$offset<-cqn.subset$glm.offset[rownames(cqn.subset$glm.offset)%in%rownames(cnts),]
  d.mont.cqn<-estimateGLMCommonDisp(d.mont,design=design)

  efit.cqn<-glmFit(d.mont.cqn,design=design)
  elrt.cqn<-glmLRT(efit.cqn,coef=2)
  tagd<-topTags(elrt.cqn,n=100,adjust="BH",p=0.05)

  #####heatmap ###########
  drawHeatmap(kexp,tags=tagd,byType="counts",cutoff=2)





  return(list(fitted=elrt.cqn,topTags=tagd,RPKM.cqn.log2=RPKM.cqn))

} #{{{ main
