#' @title produces a correlation map given wgcna network datai of specific modules and assigns a function name to the module enrichment activity
#' @description after calling wgcna the datExpr and module colors are used to create a correlation map image and a subset of these modules is plotted.
#' @param lnames this is the results from wgcna.R
#' @param read.cutoff  integer for min cutoff
#' @param plotDot boolean this plots the median correlation score for each biotypes in a line plot
#' @param recalc boolean if truen then will recalculate the bicor and pvalues
#' @param targets this is a character vector of modules of interest to subset the corMap
#' @import WGCNA
#' @import ComplexHeatmap
#' @export
#' @return images of eigengenes
wgcna_FunctionMap<-function(lnames,rnames,read.cutoff=2,recalc=FALSE,how=how,keyWord=NULL,dbName=NULL,qdbName=NULL){
  ###FIX ME: rnames has the column modules custom renamed. 
if(is.null(lnames)==TRUE){
 stop("please load wgcna.dataInput.RData")
}
###declare needed objects from load
if(is.null(lnames)==FALSE){
message(paste0("found lnames"))
} else {
 stop("did not find lnames, please run wgnca")
}
 
 df<- moduleFunctionDataFrame(db=dbName,qdb=qdbName,functionKeyWords=keyWord)
 t<-split(df,df$color)
  print(t)
 pathway<-sapply(t,function(x) sum(x$weight))
 names(pathway)<-paste0("ME",names(pathway)) ##structure must have MEcolor
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
 
  id<-match(names(pathway),rownames(moduleTraitCor)) 
  id2<-match(names(pathway),rownames(moduleTraitPvalue))
  moduleTraitCor<-moduleTraitCor[id,]
  moduleTraitPvalue<-moduleTraitPvalue[id2 ,]
 

  rTraitCor<-rnames[["traitCorRenamed"]]
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%colnames(moduleTraitCor)]
  corrMap<-bicor(t(moduleTraitCor),t(rTraitCor))
  weightedAssociation<-asinh(corrMap*(pathway^2)) #should penalize low weights and reward high weights

  ha_mix_top=HeatmapAnnotation(density_line = anno_density(weightedAssociation, 
                                            type = "line"),
                              heatmap = anno_density(weightedAssociation, 
                                            type = "heatmap"),width=unit(4,"cm"))
  weighted.activation.direction<-Heatmap(asinh(pathway),name="Activation Direction",show_row_names=FALSE)

    weight<-Heatmap(weightedAssociation, 
                  name = "asinh(weight)", 
                  top_annotation = ha_mix_top, 
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=8),
                  column_title=paste0("Weighted Module ",how))
    draw(weight+weighted.activation.direction)
   readkey()
  ##the rows are the gene Modules related to keyWord, and the columns are the module names of the repeat modules.
 #the association takes the sum of the weighted module acitivaton score for function KEYWORD, and multiplies across the correlations of the repeat modules.
  ####
###plot weighted and unweighted adjacent to each other
   pathway.unweighted<-sapply(t,function(x) sum(x$logFC))
   unweightedAssociation<-asinh(corrMap*(pathway.unweighted^2))

   ha_mix_top2=HeatmapAnnotation(density_line = anno_density(unweightedAssociation, 
                                                type = "line"),
                                 heatmap = anno_density(unweightedAssociation, 
                                                 type = "heatmap"),width=unit(3,"cm"))
  Activation.direction<-Heatmap(asinh(pathway.unweighted),name="Activation Direction",show_row_names=FALSE)
  unweight<-Heatmap(unweightedAssociation, name = "asinh(activity)", 
                   top_annotation = ha_mix_top2,
                   top_annotation_height = unit(4, "cm"),
                   column_names_gp=gpar(fontsize=8),
                   column_title=paste0(how," Activity"))

  draw(unweight+Activation.direction)
  readkey()
  dev.new(width=14,height=10)
  draw(weight+weighted.activation.direction+unweight+Activation.direction)
  readkey()

  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "");
 
 stopifnot(dim(textMatrix)==dim(moduleTraitCor))
 
 par(mar = c(6, 10, 3, 3));
 # Display the correlation values within a heatmap plot
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
  draw(weight)
  draw(unweight)
   dev.off()
 }

