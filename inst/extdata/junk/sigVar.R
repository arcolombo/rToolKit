#' @title significance variation method by Audic 
#' @description the Zhang paper used RPKM and a DE variation call using significance of variation of RPKM to rank top DE.  this method annotates the coding region and uses Audic method
#' @import biomaRt
#' @import edgeR
sigVar<-function(kexp){
 
###find RPKM
 targets<-c("Serpine1","Egln3","Scd2","Elane","Camp","Ero1l",
  "Bnip3","Pglyrp1","Ltf","Gys1","Slc2a1","lfi30","Vegfa",
  "Coro1a","Slc16a3","Prg2","Eif4ebp1","Gpi1","Arhgdib",
  "Tpi1","Mgst2","Ly6i","Cxcl5","Il6","Ccl2","lfit2",
   "Fgl2","Rsad2","Apobec1","Mmp13","Cmpk2","Ccl7","Fn1")

 counts<-collapseBundles(kexp,"gene_id")
 
 speciesMart<-useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")

 convertedEntrezID<-getBM(filters="ensembl_gene_id",
 attributes=c("ensembl_gene_id","entrezgene",speciesSymbol="mgi_symbol",
 "description"),
 values=rownames(counts),
 mart=speciesMart)
cID<-convertedEntrezID
sup<-read.csv("~/Documents/Arkas-Paper-Data/Tet2-dendritic/nature15252-s2.csv",header=TRUE,stringsAsFactors=FALSE)
colnames(sup)<-c("GeneID","log2(4h/0h)")
sup<-sup[2:nrow(sup),]

id<-convertedEntrezID$entrezgene %in%sup$GeneID

zhang<-

} ###main



