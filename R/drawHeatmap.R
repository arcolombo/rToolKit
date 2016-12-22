#' @title Draws a heatmap with annotations
#' @description draws a heatmap with annotations useful
#' @param kexp  a kallistoExperiment
#' @param tags  character of top tags from DE list
#' @param annotations boolean  to include annotation side bar, or mutation type
#' @param cutoff integer  cutoff min
#' @param tx_biotype_color_master a csv of the master color list to be reused across heatmaps for tx biotypes independent of R session
#' @param gene_biotype_color_master a csv of the master color list to be reused across heatmaps for gene biotype independent of R session.
#' @import ComplexHeatmap
#' @importFrom graphics par grid
#' @importFrom grid gpar
#' @importFrom dendsort dendsort
#' @importFrom pvclust pvclust
#' @import colorRamps
#' @export
drawHeatmap<-function(kexp,tags=NULL,annotations=TRUE,byType=c("counts","tpm"),cutoff=2,title1=NULL,tx_biotype_color_master=NULL,gene_biotype_color_master=NULL){
    stopifnot(is.null(tags)==FALSE)
    if(is.null(title1)==TRUE){
   title1="Repeat Differential Expression"
   }

   if(is.null(tx_biotype_color_master)==TRUE){
    ##generate list and print to csv for future
    tx.color.master<-repeatBiotypeColorPalette(kexp)
   write.csv(tx.color.master,file="tx_biotype_color_master.csv")
   }else{
    ##read in file, if fail to load, regenerate and print to csv
     color.df<-read.csv(tx_biotype_color_master,stringsAsFactors=FALSE,row.names=NULL,col.names=c("txb","palette"))
     ###FIX ME: have the csv enforce a certain structure.  case specific as of now.
      tx.color.master<-color.df$palette
      names(tx.color.master)<-color.df$txb
     }

 if(is.null(gene_biotype_color_master)==TRUE){
    ##generate list and print to csv for future
    gn.color.master<-repeatFamilyColorPalette(kexp)
     write.csv(gn.color.master,file="gene_biotype_color_master.csv")
   }else{
    ##read in file, if fail to load, regenerate and print to csv
     gene.color.df<-read.csv(gene_biotype_color_master,stringsAsFactors=FALSE,row.names=NULL,col.names=c("txb","palette"))
     ###FIX ME: have the csv enforce a certain structure.  case specific
      gn.color.master<-gene.color.df$palette
      names(gn.color.master)<-gene.color.df$txb
      }


   write.csv(tags,file=paste0(title1,"_drawHeatmap.repeat_top_tags.csv"))
   byType<-match.arg(byType,c("counts","tpm"))
    if(byType=="tpm"){
   rpt.tpm<-collapseTpm(kexp,"tx_id")
   } else {
   rpt.tpm<-collapseBundles(kexp,"tx_id",read.cutoff=cutoff)
   }
   df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
   rpt.targets<-rownames(tags)
  ###assign colors based on the rowRanges of kexp to be output, as a fixed table
  ##FIXME: have the annotation colors be systematic. 
 df.rpt<-df.rpt[rownames(df.rpt)%in%rpt.targets,]
  rpt.mt<-df.rpt
  rpt.mt<-as.matrix(rpt.mt)
 dev.new()
  if(nrow(rpt.mt)>=20){
  rpt.pv<-pvclust(log(1+rpt.mt),nboot=100)
  rpt.dend<-dendsort(hclust(dist(log(1+rpt.mt))),isReverse=TRUE)
  rh.rpt<-Heatmap(log(1+rpt.mt),
                  name="log(1+tpm)",
                  cluster_rows=rpt.dend,
                  cluster_columns=rpt.pv$hclust,
                  column_title=paste0(title1),
                  column_title_gp=gpar(fontsize=10),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8),
                 heatmap_legend_param=list(color_bar="continuous",
                                       #legend_direction="horizontal",
                                       legend_width=unit(5,"cm"),
                                       title_position="lefttop"))
  
 #rh.rpt2<-Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype,
 #         name="transcript_biotype",
  #        width=unit(5,"mm"))
 txb.df<-data.frame(transcript_biotype=rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype)
 txb.sel.col=list(transcript_biotype=tx.color.master)
 id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype
 rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))


 #rh.rpt3<- Heatmap(rowRanges(kexp)[rownames(rpt.mt)]$gene_biotype,
  #        name="rpt_family",
   #       width=unit(5,"mm"))

 gn.df<-data.frame(gene_biotype=rowRanges(kexp)[rownames(rpt.mt)]$gene_biotype)
 gn.sel.col=list(gene_biotype=gn.color.master)
 id3<-names(gn.sel.col[[1]])%in%rowRanges(kexp)[rownames(rpt.mt)]$gene_biotype
 rA2<-rowAnnotation(gn.df,col=list(gene_biotype=(gn.sel.col[[1]][id3])))


 draw(rh.rpt+rA+rA2)
 readkey()
 pdf(paste0(title1,"_Repeat_Heatmap_",byType,".pdf"))
  print(rh.rpt+rA+rA2)
  dev.off()
 } else {
  rh.rpt<-Heatmap(log(1+rpt.mt),
                  name="log(1+tpm)",
                 column_title=paste0("Top Repeats P.Val cut ",cutoff),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8))
   
 draw(rh.rpt)
  readkey()
 pdf(paste0(title1,"_Repeat_Heatmap_",byType,".pdf"))
  print(rh.rpt)
  dev.off()
  }

print("done.")


}#{{{main
