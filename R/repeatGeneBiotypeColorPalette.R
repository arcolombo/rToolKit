#' @title this will uniformly handle the color palette for repeat plotting
#' @description attempts to uniformly control the color assignment for biotype and family color plots for heatmapping
#' @param kexp a kallisto experiment annotated
#' @import grDevices
#' @import colorRamps
#' @import colorspace
#' @import TxDbLite
#' @import arkas
#' @export
#' @return a vector ready for rowAnnotation ComplexHeatmap
repeatGeneBiotypeColorPalette<-function(kexp){
kexp<-findRepeats(kexp)
color.df<-data.frame(txb=levels(factor(rowRanges(kexp)$gene_biotype)),palette="none",stringsAsFactors=FALSE)
  color.df[grep("DNA_element",color.df$txb),]$palette<- "#FF964D"
   color.df[grep("LINE",color.df$txb),]$palette<-"#9655E0"
     color.df[grep("LTR_element",color.df$txb),]$palette<-"#000000"
   color.df[grep("other_repeat",color.df$txb),]$palette<-"#00B734"
  color.df[grep("SINE",color.df$txb),]$palette<-"#006AC8"
  remnant.group<-length(which(color.df$palette=="none")  )
  if(remnant.group>0){
  color.df[which(color.df$palette=="none"),]$palette<-rainbow(remnant.group )
  }

 return(color.df)
} ##main
