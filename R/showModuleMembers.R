#' @title shows the genes associated with a module set
#' @description given a gene module color, we can examine the gene Ids , gene description along side an overall module enrichment activity result.  Module enrichment shows the overall enriched pathway analysis, and won't give specifics about genes.  this will display a specific module enrichment result and the module members with gene information.
#' @param lnames  the output from wgcna
#' @param biocolor a color module of interest
#' @import WGCNA
#' @export
#' @return a data frame with the eigenmodule color selected with the correlation score and pvalue and associated gene description.  downstream of this output could include filtering, and expresison plots
showModuleMembers<-function(lnames,biocolor="blue"){

if(is.null(lnames)==FALSE){
message("found wgcna.dataInput")
bwModuleColors<-lnames[["moduleColors"]]
MEs<-lnames[["MEs"]]
datExpr<-lnames[["datExpr"]]
datTraits<-lnames[["datTraits"]]
annot<-lnames[["annot"]]
} else{
stop("please run wgcna before analyszing")
} 
mm.color<-paste0("MM",biocolor)
modNames = substring(names(MEs), 3)

nGenes<-ncol(datExpr)
nSamples = nrow(datExpr);

#correlation matrix for datExpr
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
##rename for ease
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

gm.id<-which(colnames(geneModuleMembership)==mm.color)
gm<-geneModuleMembership[,gm.id]
gv.id<-which(colnames(MMPvalue)==paste0("p.MM",biocolor))
gv<-MMPvalue[,gv.id]

df<-data.frame(moduleCorrelation=gm,
               ModulePvalue=gv)

rownames(df)<-rownames(geneModuleMembership)
colnames(df)<-c(paste0("Mod.Correlation.",biocolor),
                paste0("Mod.Pvalue.",biocolor))

bm.id<-match(rownames(df),annot$ensembl_gene_id)
bm<-annot[bm.id,]
message(paste0("total NA Pvalues for ",biocolor," :",sum(is.na(df[,grep("Mod.Pvalue.",colnames(df))]))," dim: ",nrow(df) ))
df<-cbind(df,bm)
df<-df[order(df[,grep("Mod.Pvalue.",colnames(df))],decreasing=FALSE),]
return(df)

} #main
