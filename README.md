# rToolKit

"How'd you like some ice cream, Doc?"

## Work Flow
  library(arkasData)
  library(arkas)
  library(repeatToolKit)
  data(NS)
  kexp<-annotateFeatures(NS,'transcript')
  design<-model.matrix(~substring(colnames(kexp),1,1))
  rownames(design)<-colnames(kexp)
   lnames<-wgcna(kexp,whichWGCNA="single",species="Homo.sapiens",selectedPower=6,
                how="cpm",design=design,collapseBy="gene_id",annotate=TRUE)
   rnames<-wrcna(findRepeats(kexp))
   repeatTablesFromWGCNA(lnames,annotate=TRUE,byWhich="gene")
   repeatTablesFromWGCNA(rnames,annotate=FALSE,byWhich="repeat")
   qb<-qusageTablesFromWGCNA(lnames,rnames,"wgcnaDblite.sqlite","wrcnaDbLite.sqlite","qusageDbLite.sqlite")



## Exploratory

  wgcna_Heatcor(kexp)
  moduleWiseAnalysis(kexp,dbname)
  repeatWiseAnalysis(kexp)
  patientTrioPlot(kexp,printWhat="pdf")
   qusageResults(kexp,MSigDB="my.gmt")
    topHubGenes(kexp,lnames,topNumber=25)
