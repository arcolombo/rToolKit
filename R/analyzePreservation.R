#' @title analyzes the preservation of modules for a comparison colar assignment
#' @description provides a qualitative module preservation between two gene networks
#' @param geneTree1 the gene tree for the control
#' @param module1  the module color set for control
#' @param main1 descriptor
#' @param geneTree2 the gene tree for the comparison
#' @param main2 the descriptor
#' @import WGCNA
#' @export
#' @return null
analyzePreservation<-function(geneTree1=NULL,module1=NULL,main1="",geneTree2=NULL,main2="" ){

  plotDendroAndColors(geneTree1,
		     module1,
                    "Modules",
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    main=main1)
  readkey()
  plotDendroAndColors(geneTree2,
                    module1,
                    "Modules",
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    main=main2)
  readkey()
####Now to compare across kexp1 and kexp2
pdf("Compare_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTree1,
                    module1,
                    "Modules",
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    main=main1)
plotDendroAndColors(geneTree2,
                    module1,
                    "Modules",
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    main=main2)

dev.off()
  return(NULL)
} #main
