#' @title CQN and edgeR analysis
#' @description this runs CQN and edgeR on a kallistoExperiment with dynamic contrasts to select  CQN normalization across clonal stages.  This must have replicates so we are comparing samples and not specific patients.
#' @param kexp a kallistoExperiment preferably only repeats, also the repeat kexp must not be TMM Normalized
#' @param contrastMatrix a matrix with contrasts, this is specific and this function is created without fail safes
#' @param design a matrix to compare, 2 group treatment factor
#' @param comparison  the column group either pHSC , LSC ,or Blast
#' @param control  not comparison
#' @param group3   not comparison and not control
#' @import cqn
#' @import edgeR 
#' @importFrom scales alpha
#' @importFrom graphics par grid
#' @return returns 
cqnDE<-function(kexp,contrastMatrix=NULL,design=NULL,patientID=NULL, comparison,control,group3){
  stopifnot(is.null(design)==FALSE) 
  #FIX ME: currently only supports 2 group comparisons
  #FIX ME: test this out across groups Blast vs LSC,  Blast vs pHSC.   
  if(is.null(design)==TRUE){
  stopifnot(is.null(metadata(kexp)$design)==FALSE)
  design<-metadata(kexp)$design
  } 
  message("note the kexp must not be TMM normalized")
 
  dge<-DGEList(counts=counts(kexp))
  assays(kexp)$rpkm<-rpkm(dge,log=FALSE,gene.length=eff_length(kexp))
  sizeFactors.subset<-dge$samples$lib.size #grabs sequence depth
  names(sizeFactors.subset)<-rownames(dge$samples)
  uCovar<-data.frame(length=rowRanges(kexp)$tx_length,gccontent=rowRanges(kexp)$gc_content)
  rownames(uCovar)<-rownames(kexp)
  
  
  stopifnot(all(rownames(kexp)==rownames(uCovar)))
  stopifnot(all(colnames(kexp)==names(sizeFactors.subset)))

  cqn.subset<-cqn(round(assays(kexp)$rpkm),lengths=uCovar$length,x=uCovar$gccontent,sizeFactors=sizeFactors.subset,verbose=TRUE)


  par(mfrow=c(1,2))
 cqnplot(cqn.subset, n = 1, xlab = paste0("GC content Samples RPKM ",comparison," ",control," ",group3), lty = 1)
 cqnplot(cqn.subset, n = 2, xlab =paste0("length Samples RPKM ",comparison," ",control," ",group3), lty = 1)

  RPKM.cqn<-cqn.subset$y + cqn.subset$offset #on log2 scale

  RPM <- sweep(log2(counts(kexp) + 1), 2, log2(sizeFactors.subset/10^6)) #CPM
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
 ylim = c(-4,4), xlim = c(0,19), col = "red")
 points(A.cqn[whLow], M.cqn[whLow], cex = 0.55, pch = 16, col = "blue")


########edgeR differential expression analysis

d.mont <- DGEList(counts = counts(kexp), lib.size = sizeFactors.subset, group = sapply(strsplit(pData(kexp)$ID,"_"),function(x) x[1])
, genes = uCovar)

  ##FIX ME add multi group
  d.mont$offset<-cqn.subset$glm.offset
  d.mont.cqn<-estimateGLMCommonDisp(d.mont,design=design)

  ###FIX ME::: have users input this
 # if(is.null(contrastMatrix)==TRUE){
 # contrast.matrix<-makeContrasts(LSCvpHSC=LSC-pHSC,Blast-LSC,Blast-pHSC,pHSC-LSC-Blast,levels=groupModel)
 # } else {
 #  contrast.matrix<-contrastMatrix
 #  print(contrast.matrix)
# }  ###issue with dispersion estimation with multi groups,  getting LAPACK errors... 2 group only for now....
  
  efit.cqn<-glmFit(d.mont.cqn,design=design)
  elrt.cqn<-glmLRT(efit.cqn,coef=2)
  tagd<-topTags(elrt.cqn,n=100,adjust="BH",p=0.05)
  return(list(fitted=elrt.cqn,topTags=tagd))

} #{{{ main
