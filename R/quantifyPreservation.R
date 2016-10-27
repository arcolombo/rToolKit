#' @title quantifies preservation of modules
#' @description returns a Z score summary of the module preservation. Performs permuatation testing of the modules. This measures how well the preservation of grouping of modules of the control match in the comparison.
#' @param datExpr1 is the control dat expression
#' @param datExpr2 is the comparison expression, the network that we will test for conservation
#' @param modules is the color module vector in the control
#' @import WGCNA
#' @export
#' @return z scores
quantifyPreservation<-function(datExpr1=NULL,datExpr2=NULL,modules=NULL){


  multiExpr = list(A1=list(data=t(datExpr1)),A2=list(data=t(datExpr2)))
  multiColor = list(A1 = modules)
  mp=modulePreservation(multiExpr,
                        multiColor,
                        referenceNetworks=1,
                        verbose=3,
                        networkType="signed",
                        nPermutations=30,
                        maxGoldModuleSize=100,
                        maxModuleSize=400)
  stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
  return(stats[order(-stats[,2]),c(1:2)] )

  
} #main
