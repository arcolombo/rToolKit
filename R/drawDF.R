#' @title draws a heatmap downstream of PoissonSeq
#' @description draws a heatmap from the poissonseq DE list
#' @param kexp a kallisto experiment usually a repeat stage level
#' @param res  poissonSeq fit
#' @param p.value numeric floor 
#' @import ComplexHeatmap
#' @import pvclust
#' @import dendsort
#' @export
#' @return images and stuffs
drawDF<-function(kexp,res=res,p.value=0.05,cutoff=1){

rpt.targets<-res$gname[which(res$fdr<=p.value)]
###TPM Heatmap of top fitted
df.rpt<-collapseTpm(kexp,"tx_id")
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


}###main
