#' @title Draws a heatmap with annotations
#' @description draws a heatmap with annotations useful
#' @param kexp  a kallistoExperiment
#' @param tags  character of top tags from DE list
#' @param annotations boolean  to include annotation side bar, or mutation type
#' @param cutoff integer  cutoff min
#' @import ComplexHeatmap
#' @importFrom graphics par grid
#' @importFrom grid gpar
#' @importFrom dendsort dendsort
#' @importFrom pvclust pvclust
#' @export
drawHeatmap<-function(kexp,tags=NULL,annotations=TRUE,byType=c("counts","tpm"),cutoff=2){
    stopifnot(is.null(tags)==FALSE)
   byType<-match.arg(byType,c("counts","tpm"))
    if(byType=="tpm"){
   rpt.tpm<-collapseTpm(kexp,"tx_id")
   } else {
   rpt.tpm<-collapseBundles(kexp,"tx_id",read.cutoff=cutoff)
   }
   df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
   rpt.targets<-rownames(tags)

 df.rpt<-df.rpt[rownames(df.rpt)%in%rpt.targets,]
  rpt.mt<-df.rpt
  rpt.mt<-as.matrix(rpt.mt)

  if(nrow(rpt.mt)>=20){
  rpt.pv<-pvclust(log(1+rpt.mt),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+rpt.mt))),isReverse=TRUE)
  rh.rpt<-Heatmap(log(1+rpt.mt),
                  name="log(1+tpm)",
                  cluster_rows=rpt.dend,
                  cluster_columns=rpt.pv$hclust,
                  column_title=paste0("Top Repeats P.Val cut ",cutoff),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8))
 rh.rpt2<-Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm"))
 rh.rpt3<- Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
 draw(rh.rpt+rh.rpt2+rh.rpt3)
 } else {
  rh.rpt<-Heatmap(log(1+rpt.mt),
                  name="log(1+tpm)",
                 column_title=paste0("Top Repeats P.Val cut ",cutoff),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8))
   }
readkey()
}#{{{main