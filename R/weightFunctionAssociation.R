#' @title given the module of interest printed to a csv, we can read that data, and analyze the pathway behavior of specific pathways
#' @description reads the qusage module pathway activity and analyzes specific behvaviors within that function
#' @param db the wgcnaDbName
#' @param qdb the qusageDbLite name
#' @param functionKeyWords this is a single character that will be grep'd from the list of pathway modules 
#' @export
#' @return a data frame with the queried keyword in every module and pathway information
weightFunctionAssociations<-function(lnames,rnames,read.cutoff=2,recalc=FALSE,how=how,dbName=NULL,qdbName=NULL,keyWords=NULL){
 
   df<-pickPathway(qusageDbLite(qdb),keyWord=keyWords)
  print(df)
  ###This requires that rnames has the columns renamed txBiotypes that lead the module
  ### this ranks the correlation matrix by the rank of the pathway function queried.  this rewards correlations to highly ranked pathway query functions, and should penalize slighly low ranked correlations to low ranked pathway query functions.  normalizes the module pathway size.
  ### averages the logFC in the query, assumes that the query will return multiple similar pathways, so averages the queries into an average function logFC (activity direction) 
  ranking.weight<-sapply(df,function(x) x[grep("Total",rownames(x)),]$ranking/x[grep("Total",rownames(x)),]$module.size)
  activation.direction<-sapply(df,function(x) x[grep("Total",rownames(x)),]$logFC/x[grep("Total",rownames(x)),]$query.size)


if(is.null(lnames)==TRUE){
 stop("please load wgcna.dataInput.RData")
}
###declare needed objects from load
if(is.null(lnames)==FALSE){
message(paste0("found lnames"))
} else {
 stop("did not find lnames, please run wgnca")
}



  names(ranking.weight)<-paste0("ME",names(ranking.weight)) ##structure must have MEcolor
  bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
# Will display correlations and their p-values
# Define numbers of genes and samples
  nGenes= ncol(datExpr);
  nGenes= ncol(datExpr);
  nSamples = nrow(datExpr);
 # Recalculate MEs with color labels
  if(recalc==TRUE){
   moduleTraitCor = bicor(MEs, datTraits, use = "all.obs");
   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
   } else{
   moduleTraitCor<-lnames[["moduleTraitCor"]]
   moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
   }

 id<-match(names(ranking.weight),rownames(moduleTraitCor))
  id2<-match(names(ranking.weight),rownames(moduleTraitPvalue))
  moduleTraitCor<-moduleTraitCor[id,]
  moduleTraitPvalue<-moduleTraitPvalue[id2 ,]


  rTraitCor<-rnames[["traitCorRenamed"]] 
   if(ncol(as.data.frame(moduleTraitCor))>1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%colnames(moduleTraitCor)]
  }else if(ncol(as.data.frame(moduleTraitCor))==1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%names(moduleTraitCor)]
  }
  corrMap<-bicor(t(moduleTraitCor),t(rTraitCor))
  weightedAssociation<-(corrMap*(ranking.weight)) #should penalize low weights and reward high weights

  ha_mix_top=HeatmapAnnotation(density_line = anno_density(weightedAssociation,
                                            type = "line"),
                              heatmap = anno_density(weightedAssociation,
                                            type = "heatmap"),width=unit(4,"cm"))
  weighted.activation.direction<-Heatmap(asinh(activation.direction*ranking.weight),name="Activation Direction",show_row_names=FALSE)

    weight<-Heatmap(asinh(weightedAssociation),
                  name = "asinh(weight)",
                  top_annotation = ha_mix_top,
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=8),
                  column_title=paste0("Weighted Module ",how))
    draw(weight+weighted.activation.direction)
   readkey()
 ###SO FAR SO GOOD>>> Fix rest.

 
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "");

 stopifnot(dim(textMatrix)==dim(moduleTraitCor))

 par(mar = c(6, 10, 3, 3));
 plot.new()
  labeledHeatmap(Matrix=moduleTraitCor,
,                xLabels=names(datTraits),
                 yLabels= rownames(moduleTraitCor),
                 ySymbols=rownames(moduleTraitCor),
                 colorLabels=FALSE,
               colors=greenWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               yColorWidth=0.07,
               cex.lab.y=.7,
               colors.lab.y=1.3,
               main = paste0("Module-Repeat ",how," Biotype relationships"))


  readkey()


   pdf(paste0("weighted_correlation_",how,"_plots.pdf"),width=15,height=10)
    par(mar = c(6, 10, 3, 3));
    labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.3,
                zlim = c(-1,1),
                yColorWidth=0.07,
                cex.lab.y=.7,
                colors.lab.y=1.3,
                main = paste0("Module-Repeat ",how," Biotype relationships"))
      draw(weight+weighted.activation.direction)
   dev.off()
 
}
