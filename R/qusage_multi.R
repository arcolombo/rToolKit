#' @title run qusage over multiple gene sets
#' @description since the WGCNA enrichment algorithm is 'experimental' it can only provide preliminary analyses; additionally it returns warnings regarding pvalue calculations on correlations.  qusage can be run to verify using KEGG, GO, and immuneSignatures.  Note that WGCNA:Go enrichment inputs Entrez IDs which is gene specific, so in this function we collapse by gene ids, and run qusage 
qusage_multi<-function(kexp,labels,dbDir="~/Documents/Arkas-Paper-Data/MSigDB/",db="c2.cp.kegg.v5.1.symbols.gmt",comparison="LSC",control="pHSC"){


kexp<-kexp2Group(kexp,comparison=comparison,control=control)

eset<-collapseBundles(kexp,"gene_name")
message("creating labels")
compID<-grep(comparison,colnames(eset))
controlID<-grep(control,colnames(eset))
labels<-colnames(eset)[c(compID,controlID)]
##cut the patient ID off
labels<-sapply(sapply(labels,function(x) strsplit(x,"_")),function(x) x[1])
labels<-as.vector(labels)
contrast<-paste0(comparison,"-",control)
#FIX ME: generalize this method to run over all DB.
MSIG.geneSets<-read.gmt(paste0(dbDir,db))
qs.results.msig = qusage(eset, labels, contrast, MSIG.geneSets)
print(numPathways(qs.results.msig))
p.vals = pdf.pVal(qs.results.msig)



} #Main
