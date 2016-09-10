#' @title analyzes kexp data
#' @description analyzes kexp repeat and gene data but requires a design matrix.  uses RUV and a limma/voom model linear
#' @import arkas
#' @import RUVSeq
#' @import ComplexHeatmap
#' @param kexp  mpn vs Mx Cre lines
#' @export
kexpAnalysis<-function(kexp,byWhich=c("gene","repeat"),read.cutoff=1,adjustBy=c("BH","none")){
#kexp has design matrix and is annotated
#running GWA and RWA
   adjustBy<-match.arg(adjustBy,c("BH","none"))
   design<-metadata(kexp)$design
   if(byWhich=="gene"){
   gwa<-geneWiseAnalysis(kexp,design=metadata(kexp)$design,how="cpm",read.cutoff=read.cutoff,species="Mus.musculus",adjustBy="BH")
  
  ruvGWA<-ruvNormalization(kexp,byLevel="gene_id")
  design<-cbind(design,ruvGWA$RUV$W)
  norm.gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",read.cutoff=read.cutoff,species="Mus.musculus",adjustBy="BH")

  tpm<-collapseTpm(kexp,"gene_id")
   rpt.mt<-tpm[rownames(tpm)%in%norm.gwa$topGenes,]
   rh<-Heatmap(rpt.mt,
               column_names_gp=gpar(fontsize=7),
               row_names_gp=gpar(fontsize=10),
                name="tpm",
                column_title="RUV TPM")
  rh2<-Heatmap(rowRanges(kexp)[norm.gwa$topGenes]$tx_biotype,name="Tx_biotype")
  rh3<-Heatmap(rowRanges(kexp)[norm.gwa$topGenes]$gene_biotype,name="Gene_biotype")
  draw(rh+rh2+rh3)
  readkey()

  
  return(norm.gwa)
  } else {

   rwa<-repeatWiseAnalysis(kexp,
                          design=design,
                          how="cpm",
                          species="Mus.musculus",
                          adjustBy=adjustBy,
                          read.cutoff=read.cutoff)
   tpm<-collapseTpm(kexp,"tx_id")
   rpt.mt<-tpm[rownames(tpm)%in%rwa$topGenes,]
   rh<-Heatmap(asinh(rpt.mt),
               column_names_gp=gpar(fontsize=7),
               row_names_gp=gpar(fontsize=10),
                name="asinh(tpm)",
                column_title="RUV TPM")
  rh2<-Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype,name="Tx_biotype",width=unit(5,"mm"))
  rh3<-Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$gene_biotype,name="Gene_biotype",width=unit(5,"mm"))
  draw(rh+rh2+rh3)
  readkey()
  
  ###FIX ME:: add ensmblToHUGO call 
 
  #look for Cre (should be up) DNMT3A should be 50% down wrt MPN  erythocytes, lymphocytes and lymphoids should be down in MxCre
   return(rwa)
  } #else repeats

} #main

