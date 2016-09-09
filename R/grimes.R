#' @title analyzes grimes data
#' @description analyzes grimes repeat data
#' @import arkas
#' @import RUVSeq
#' @import ComplexHeatmap
#' @param grimes  mpn vs Mx Cre lines
grimes<-function(grimes){
#grimes has design matrix and is annotated
#running GWA and RWA
   design<-metadata(grimes)$design
   gwa<-geneWiseAnalysis(grimes,design=metadata(grimes)$design,how="cpm",read.cutoff=1,species="Mus.musculus",adjustBy="BH")
  
  ruvGWA<-ruvNormalization(grimes,byLevel="gene_id")
  design<-cbind(design,ruvGWA$RUV$W)
  norm.gwa<-geneWiseAnalysis(grimes,design=design,how="cpm",read.cutoff=2,species="Mus.musculus",adjustBy="BH")
  rwa<-repeatWiseAnalysis(grimes,design=design,how="cpm",species="Mus.musculus",adjustBy="BH",read.cutoff=2)
   tpm<-collapseTpm(grimes,"tx_id")
   rpt.mt<-tpm[rownames(tpm)%in%rwa$topGenes,]
   rh<-Heatmap(rpt.mt,
               column_names_gp=gpar(fontsize=7),
               row_names_gp=gpar(fontsize=10),
                name="tpm",
                column_title="RUV TPM")
  rh2<-Heatmap(rowRanges(grimes)[rwa$topGenes]$tx_biotype,name="Tx_biotype")
  rh3<-Heatmap(rowRanges(grimes)[rwa$topGenes]$gene_biotype,name="Gene_biotype")
  draw(rh+rh2+rh3)

}
