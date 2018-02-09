# rToolKit

"How'd you like some ice cream, Doc?"

## Work Flow

   lnames<-wgcna(kexp)
   rnames<-wrcna(findRepeats(kexp))
   repeatTablesFromWGCNA(lnames,annotate=TRUE,byWhich="gene")
   repeatTablesFromWGCNA(rnames,annotate=FALSE,byWhich="repeat")
   qb<-qusageTablesFromWGCNA(lnames,rnames,"wgcnaDblite.sqlite","wrcnaDbLite.sqlite","qusageDbLite.sqlite")



## Exploratory

  wgcna_Heatcor(kexp)
  moduleWiseAnalysis(kexp,dbname....)
  repeatWiseAnalysis(kexp..)
  patientTrioPlot(kexp,printWhat="pdf")
  drawBoxPlots(kexpByTrio(kexp),numberComparisons=2,read.cutoff=1,adjustBy="none",wilcox.Alternative="greater",testMedians=TRUE,title1="Absolute Differential Transposable Element Activity Comparing ",ylab1=ylab1,xlab1=xlab1,xlab2=xlab2)
   qusageResults(kexp,MSigDB="my.gmt",...)
   topHubGenes(kexp,lnames,topNumber=25)
