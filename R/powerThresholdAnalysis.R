#' @title this will analysze the power threshold beta
#' @description although the network is robust to the choice of beta, this method will call the pickSoftThreshold values and print to a pdf for future use.
#' @param lnames this is the list from wgcna.R
#' @param datExpr this is the log2 XR rows as samples columns as genes/repeats
#' @param how a character param for pdf print
#' @import WGCNA
#' @export
powerThresholdAnalysis<-function(lnames=NULL,datExpr=NULL,how="gene"){
 
 if(is.null(datExpr)==TRUE){
  stopifnot(is.null(lnames)==FALSE)
 datExpr<-lnames[["datExpr"]]
 }else if(is.null(lnames)==TRUE){
  stopifnot(is.null(datExpr)==FALSE)
  datExpr<-datExpr
 }else{
 stop("Please enter lnames which is the output of wgcna.R, or t(eset) not both...\n")
  }
 powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  readkey()
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
   text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers,cex=cex1,col="red");
  selectedPower<-readPower()
##FIX ME::: THIS IS NOT PRINTING PDF PROPERLY ? COMES OUT BLANK
 pdf(paste0("WGCNAselectedPower_",how,"analysis.pdf"),width=11,height=9)
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
   plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
   text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers,cex=cex1,col="red");
  dev.off()

return(selectedPower)

} ##main
