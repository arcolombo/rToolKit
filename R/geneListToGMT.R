#' @title convert list of genes to GMT format for qusage calls
#' @description given a gene list write to a .gmt file
#' @param geneList a list of genes with names(list) 
#' @param fileName output.gmt
#' @param citation the source of gene list
#' @export
geneListToGMT<-function(geneList=NULL,fileName="test.gmt",citation="www.na.com"){
   for(i in 1:length(names(geneList))){
    x<-data.frame(c(citation,geneList[[i]]))
    colnames(x)<-names(geneList)[i]
   ##there is a bugger where the last gene prints an eol as tab and \n is needed for the last gene
    xnrow<-nrow(x)
    x.first<-seq(1,(xnrow-1))

   write.table(x[x.first,],
               file=fileName,
               append= T, 
               sep='\t',
               quote=FALSE,
               row.names=FALSE,
               col.names=colnames(x),
               eol="\t" )
   write.table(x[xnrow,],
               file=fileName,
               append=T,
               sep='\t',
               quote=FALSE,
               row.names=FALSE,
               col.names=FALSE,
               eol="\n")
  }#loop
}
