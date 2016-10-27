#' @title this plots all the dendrograms from a net or bwnet call
#' @description for wgcna.R call, the network dendrograms are clusterd based on average expression, one for each module. this iteratively goes through each and plots them.
#' @param bwnet$dendrograms this is a gene tree from the block net call
#' @param net$dendrograms  a gene tree from single net call
#' @param whichWGCNA single or block level analysis, block is a bit less rigid
#' @export
#' @return cool network images and stuff
wgcna_plotAll_dendrograms<-function(bwnet=NULL,whichWGCNA=c("single","block"),bwModuleColors=NULL,bwLabels=NULL,how="cpm",byWhich="gene"){
  if(is.null(bwLabels)==TRUE){
  bwLabels<-bwnet$colors
  }
  bwLabels = matchLabels(bwnet$colors,bwnet$colors)
  if(is.null(bwModuleColors)==TRUE){
  bwModuleColors = labels2colors(bwLabels) 
  }
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  if(whichWGCNA=="block") {
   for(i in 1:length(bwnet$dendrograms)){
     plotDendroAndColors(bwnet$dendrograms[[i]],
                      bwModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors",
                     main =paste0(byWhich," dendrogram ",how," module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
   readkey()
    }
   pdf(paste0("allBlockDendroGrams_",byWhich,"_",how,".pdf"),width=12,height=12)
     for(i in 1:length(bwnet$dendrograms)){
     plotDendroAndColors(bwnet$dendrograms[[i]],
                      bwModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors",
                     main =paste0(byWhich," dendrogram ",how," module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
     }
  dev.off()
  } else {

 for(i in 1:length(bwnet$dendrograms)){
     plotDendroAndColors(bwnet$dendrograms[[i]],
                      bwModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors",
                     main =paste0(byWhich," dendrogram ",how," and module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
   readkey()
    }

  pdf(paste0("allSingleDendroGrams_",byWhich,"_",how,".pdf"),width=12,height=12)
   for(i in 1:length(bwnet$dendrograms)){
     plotDendroAndColors(bwnet$dendrograms[[i]],
                      bwModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors",
                     main =paste0(byWhich," dendrogram ",how," and module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
     }
  dev.off()
  }

} #main
