# rToolKit
This is an rToolKit
#Work Flow
lnames<-wgcna(kexp)
rnames<-wrcna(findRepeats(kexp))
repeatTablesFromWGCNA(lnames,annotate=TRUE,which="gene")
repeatTablesFromWGCNA(rnames,annotate=FALSE,which="repeat")
qb<-qusageTablesFromWGCNA(lnames,rnames,"wgcnaDblite.sqlite","wrcnaDbLite.sqlite","qusageDbLite.sqlite")

##Exploratory
wgcna_Heatcor(kexp)
moduleWiseAnalysis(kexp,dbname....)
repeatWiseAnalysis(kexp..)
