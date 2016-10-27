#' @title this analyzes tet2 and dnmt3a dKO HSC
#' @description title
#' @param kexp  this is the full kexp from peggy's HSC, we will split into each group and run a linear model regression on the data for repeats
peggyTet2Analysis<-function(kexp){

design<-metadata(kexp)$design
##make 2 designs
tetid<-grep("tet2KO",colnames(kexp))
hsc<-grep("wthsc",colnames(kexp))
tet2<-kexp[,c(tetid,hsc)]
aID<-grep("3aDKO",colnames(kexp))
aDKO<-kexp[,c(aID,hsc)]

 tet_design<-design[rownames(design)%in%colnames(tet2),]
 aDKO_design<-design[rownames(design)%in%colnames(aDKO),]
  metadata(tet2)$design<-tet_design
  metadata(aDKO)$design<-aDKO_design
  tet_rwa<-kexpAnalysis(tet2,byWhich="repeat",adjustBy="none")
    aDKO_rwa<-kexpAnalysis(aDKO,byWhich="repeat",adjustBy="none")
  peggy<-kexpAnalysis(kexp,byWhich="repeat",adjustBy="none")

}
