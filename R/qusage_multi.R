#' @title run qusage over multiple gene sets
#' @description since the WGCNA enrichment algorithm is 'experimental' it can only provide preliminary analyses; additionally it returns warnings regarding pvalue calculations on correlations.  qusage can be run to verify using KEGG, GO, and immuneSignatures.  Note that WGCNA:Go enrichment inputs Entrez IDs which is gene specific, so in this function we collapse by gene ids, and run qusage 
#' @param kexp a 2 stage level kexp
#' @param dbDir path to MSIGDB
#' @param db  char name of db
#' @param comparison character of Numerator
#' @param control  denominator
#' @param points.max  number of values to distribute into tdist
#' @param table.number number to report
#' @import qusage
#' @export
#' @return qusage table 
qusage_multi<-function(kexp,dbDir="~/Documents/Arkas-Paper-Data/MSigDB/",db="c2.cp.kegg.v5.1.symbols.gmt",comparison="LSC",control="pHSC",points.max=2^12,table.number=15){

####FIX ME : add patient specific effects
kexp<-kexp2Group(kexp,comparison=comparison,control=control)

eset<-collapseBundles(kexp,"gene_name")
eset.log2<-log2(1+eset)
message("creating labels")
compID<-grep(comparison,colnames(eset.log2))
controlID<-grep(control,colnames(eset.log2))
labels<-colnames(eset.log2)[c(compID,controlID)]
##cut the patient ID off
labels<-sapply(sapply(labels,function(x) strsplit(x,"_")),function(x) x[1])
labels<-as.vector(labels)
contrast<-paste0(comparison,"-",control)
#FIX ME: generalize this method to run over all DB.
MSIG.geneSets<-read.gmt(paste0(dbDir,db))
qs.results.msig = qusage(eset.log2, labels, contrast, MSIG.geneSets,n.points=points.max)
print(numPathways(qs.results.msig))
p.vals = pdf.pVal(qs.results.msig)
q.vals = p.adjust(p.vals, method="fdr")
qsTable(qs.results.msig, number=table.number)
readkey()
plot(qs.results.msig)
readkey()
return(qs.results.msig)
} #Main
