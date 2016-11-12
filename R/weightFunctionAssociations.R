#' @title given the module of interest printed to a csv, we can read that data, and analyze the pathway behavior of specific pathways
#' @description reads the qusage module pathway activity and analyzes specific behvaviors within that function
#' @param db the wgcnaDbName
#' @param qdb the qusageDbLite name
#' @param functionKeyWords this is a single character that will be grep'd from the list of pathway modules 
#' @import ComplexHeatmap
#' @import WGCNA
#' @import circlize
#' @export
#' @return a data frame with the queried keyword in every module and pathway information
weightFunctionAssociations<-function(lnames,rnames,recalc=FALSE,how=how,dbName=NULL,qdbName=NULL,keyWords=NULL,write.out=FALSE){
 
   df<-pickPathway(qusageDbLite(qdbName),keyWord=keyWords)
   csvName<-paste0(gsub(" ","_",how),"_functionAssoc.csv")
  if(write.out==TRUE){
   lapply(df,function(x) write.table(data.frame(x),file=csvName,append=TRUE,sep=',',row.names=FALSE,col.names=TRUE))
   }
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
   names(activation.direction)<-paste0("ME",names(activation.direction)) ##structure must have MEcolor

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
  moduleTraitCor<-moduleTraitCor[!is.na(rownames(moduleTraitCor)),]
  moduleTraitPvalue<-moduleTraitPvalue[!is.na(rownames(moduleTraitPvalue)),]

  rTraitCor<-rnames[["traitCorRenamed"]] 
   if(ncol(as.data.frame(moduleTraitCor))>1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%colnames(moduleTraitCor)]
  }else if(ncol(as.data.frame(moduleTraitCor))==1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%names(moduleTraitCor)]
  }
  corrMap<-bicor(t(moduleTraitCor),t(rTraitCor))
  corrMap.pvalue<-corPvalueStudent(corrMap,nrow(corrMap))
  ranking.weight<-ranking.weight[names(ranking.weight)%in%rownames(corrMap)]
  activation.direction<-activation.direction[names(activation.direction)%in%rownames(corrMap)]
  weightedAssociation<-(corrMap*(ranking.weight)) #should penalize low weights and reward high weights
 ###weighted correlations 
   ###
 #######this uses the correlation test statistic however the weighted correlations are not [0,1] where the positive weights reward correlations.  so we take the abs(1-r^2/var(rank)) which does not return complex numbers.
 ##the problem with abs(1-r^2) is that for values >1 it levels off at 2, and may miss these weighted values.
 ##so by taking abs(1-weight^2/var(ranking.weight))  it pushes the asymptote out by the sd(ranking.weight), thus if the weights are noisy the asympotote shifts further out. which is more conservative.
 ##if the keyWord query is small, then sd(rank.weight) will be smaller leading to higher FPR.
 ###this means that for values weighted closer to the sd(ranking.weight) these will have higher p.values so the higher pvalues will be those positively effected by the weights, and it will penalize according to the weights.  where before any r value close to 1 lands in critical region.  weighted pvalues will land in critical regions for any value close to the sd(ranking.weight) which will tend to be value rewarded by weights.
  
  print(paste0("weighted center:",sd(ranking.weight)))
   asymp.T<-sqrt(nrow(weightedAssociation)-2)*weightedAssociation/(sqrt(abs(1-weightedAssociation^2/var(ranking.weight)))) 
   weighted.pvalue<- 2 * pt(abs(asymp.T), nSamples - 2, lower.tail = FALSE)

   pdf(paste0("Method_",how,"_Comparisons_Weighted_Correlations.pdf"))
   asymp.T<-sqrt(nrow(weightedAssociation)-2)*weightedAssociation/(sqrt(abs(1-weightedAssociation^2/var(ranking.weight))))
   weighted.pvalue<- 2 * pt(abs(asymp.T), nSamples - 2, lower.tail = FALSE)
   hist(weighted.pvalue,main="weighted pvalues")
 
  large<-which(weightedAssociation>1)
  small<-which(weightedAssociation< ( -1) )
  wind<-weightedAssociation
  wind[large]<-1
  wind[small]<-(-1)  
  weighted.pvalue3<-corPvalueStudent(wind,nrow(wind))
  hist(weighted.pvalue3,main="windsorized weighted pvalues")
  hist(corrMap.pvalue,main="Default Standard (Unweighted) |p|<=1")
  dev.off()
 ####################
  ###winsor method... does okay
  # large<-which(weightedAssociation>1)
  # small<-which(weightedAssociation< ( -1) )
  # wind<-weightedAssociation
  # wind[large]<-1
  # wind[small]<-(-1)  
  # weighted.pvalue<-corPvalueStudent(wind,nrow(weightedAssociation))
#############

  ha_mix_top=HeatmapAnnotation(density_line = anno_density(weightedAssociation,
                                            type = "line"),
                              heatmap = anno_density(weightedAssociation,
                                            type = "heatmap"),width=unit(4,"cm"))
  weighted.activation.direction<-Heatmap(asinh(activation.direction*ranking.weight),name="Activation Direction",show_row_names=FALSE)

    weight<-Heatmap(asinh(weightedAssociation),
                   col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  name = "asinh(weight)",
                  top_annotation = ha_mix_top,
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=8),
                  column_title=paste0("Weighted Module ",how),
                  cell_fun=function(j,i,x,y,w,h,col){
                  asymp.T<-sqrt(nrow(weightedAssociation)-2)*weightedAssociation/sqrt(abs(1-weightedAssociation^2))
                 weighted.pvalue<- 2 * pt(abs(asymp.T), nSamples - 2, lower.tail = FALSE)
                  if(weighted.pvalue[i,j]<0.06){
                   grid.text(sprintf("%.3f", weighted.pvalue[i,j]),x,y)
                  }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
                  })


  
  draw(weight+weighted.activation.direction)
 readkey()
 ###########
 ###unweighted associations 

 unweighted_mix_top=HeatmapAnnotation(density_line = anno_density(corrMap,
                                            type = "line"),
                              heatmap = anno_density(corrMap,
                                            type = "heatmap"),width=unit(4,"cm"))

 unweight<-Heatmap(asinh(corrMap),
                  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  name = "asinh(cor)",
                  top_annotation = unweighted_mix_top,
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=8),
                  column_title=paste0("unWeighted Module ",how),
                 cell_fun=function(j,i,x,y,w,h,col){
                  if(corPvalueStudent(corrMap,nrow(corrMap))[i,j]<0.06){
                   grid.text(sprintf("%.3f", corPvalueStudent(corrMap,nrow(corrMap))[i,j]),x,y)
                  }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
                  })
   unweighted.activation.direction<-Heatmap(asinh(activation.direction),name="Activation Direction",show_row_names=FALSE)



    draw(unweight+unweighted.activation.direction)
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
    draw(unweight+unweighted.activation.direction)
   dev.off()
 
}
