#' @title given the module of interest printed to a csv, we can read that data, and analyze the pathway behavior of specific pathways
#' @description reads the qusage module pathway activity and analyzes specific behvaviors within that function
#' @param db the wgcnaDbName
#' @param qdb the qusageDbLite name
#' @param functionKeyWords this is a single character that will be grep'd from the list of pathway modules 
#' @import ComplexHeatmap
#' @import WGCNA
#' @import circlize
#' @import gridExtra
#' @import ggplot2
#' @import fitdistrplus
#' @import logspline
#' @export
#' @return a data frame with the queried keyword in every module and pathway information
weightFunctionAssociations<-function(lnames,rnames,recalc=FALSE,how=how,dbname=NULL,qdbname=NULL,keyWords=NULL,write.out=FALSE){
    dbName<-dbname
    qdbName<-qdbname
   df<-pickPathway(qusageDbLite(qdbname),keyWord=keyWords)
   csvName<-paste0(gsub(" ","_",how),"_functionAssoc.csv")
  if(write.out==TRUE){
   lapply(df,function(x) write.table(data.frame(x),file=csvName,append=TRUE,sep=',',row.names=FALSE,col.names=TRUE))
   }
  print(df)
  ###This requires that rnames has the columns renamed txBiotypes that lead the module
  ### this ranks the correlation matrix by the rank of the pathway function queried.  this rewards correlations to highly ranked pathway query functions, and should penalize slighly low ranked correlations to low ranked pathway query functions.  normalizes the module pathway size.
 
  ranking.weight<-sapply(df,function(x) x[grep("Total",rownames(x)),]$ranking/x[grep("Total",rownames(x)),]$module.size)
  activation.direction<-sapply(df,function(x) as.numeric(x[grep("Total",rownames(x)),]$logFC)/as.numeric(x[grep("Total",rownames(x)),]$query.size))


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
   if(is.null(rnames[["traitCorRenamed"]])==FALSE){
  rTraitCor<-rnames[["traitCorRenamed"]] 
  }else{
  rTraitCor<-rnames[["moduleTraitCor"]]
  }
   if(ncol(as.data.frame(moduleTraitCor))>1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%colnames(moduleTraitCor)]
  }else if(ncol(as.data.frame(moduleTraitCor))==1){
  rTraitCor<-rTraitCor[,colnames(rTraitCor)%in%names(moduleTraitCor)]
  }
  corrMap<-bicor(t(moduleTraitCor),t(rTraitCor))
  corrMap.pvalue<-corPvalueStudent(corrMap,ncol(corrMap))
  corrMap.pvalue[which(is.na(corrMap.pvalue))]<-1
  ranking.weight<-ranking.weight[names(ranking.weight)%in%rownames(corrMap)]
  activation.direction<-activation.direction[names(activation.direction)%in%rownames(corrMap)]
 ##instead of using the ranking.weight raw score find out where they fit onto a normal distribution. and scale appropriately
  ##the log(data) helps normalize the data transform into normality.  
  ##FIX ME: ensure/enforce that the data transform works? 
  
 normal.scaled.weights<-pnorm(log(ranking.weight),mean=mean(log(ranking.weight)),sd=sd(log(ranking.weight)),lower.tail=TRUE)
  #print(shapiro.test(normal.scaled.weights))
  ####optimize the weights to be normal and the weightedAssoc to be normal
  ###loop through the bottom, top pairs of percentils, testing normality, while adjusting both ends balanced

 adj.per<-normal.scaled.weights
 for(i in 1:as.integer(length(order(rank(adj.per)))/2)){

   max.id<-order(rank(adj.per),decreasing=TRUE)[i]
   min.id<-order(rank(adj.per),decreasing=FALSE)[i]
   message(paste0("i ",i))
   for(k in 2:6){
    X<-adj.per
   message(paste0("k ",k))
   if(i==1){
    ##100% and 0% quantile
    step1<-abs(X[max.id ]-1)/k
    step2<-abs(X[min.id ]-0)/k
   }else if(i>1){ 
    ##middle percentiles
    step1<-abs(X[max.id ]- X[order(rank(X),decreasing=TRUE)[i-1]] )/k
    step2<-abs(X[min.id ]- X[order(rank(X),decreasing=FALSE)[i-1]] )/k 
       print(paste0("step2: ",step2))
    }

 X[max.id]<-X[max.id]+step1
  X[min.id]<-abs(X[min.id]-step2)
  print(X[max.id])
  print(X[min.id])
  if(shapiro.test(X)$p.value<0.06 && shapiro.test(corrMap*X)$p.value<0.06) {
   message(paste0(names(X)[i]," at ",k," weight was not normal retrying..."))
   next #next k iteration
   }else if(shapiro.test(X)$p.value>0.06 && shapiro.test(corrMap*X)$p.value>0.05){
  message("corrMap*adjusted weights are normal")
    adj.per[max.id]<-X[max.id]
  adj.per[min.id]<-X[min.id]
   print(c(adj.per[max.id],adj.per[min.id]))
    break
     }
   }##loop k partition
 } #quantile loop
 normal.scaled.weights<-adj.per 

 weightedAssociation<-(corrMap*(normal.scaled.weights)) #should penalize low weights and reward high weights
  if(write.out==TRUE){
  write.csv(weightedAssociation,file=paste0("weightedCorrelations_",gsub(" ","_",how),".csv"))
  write.csv(normal.scaled.weights,file=paste0("weights_",gsub(" ","_",how),".csv"))
  write.csv(ranking.weight,file=paste0("rankingActivity_",gsub(" ","_",how),".csv"))
  }
 ###weighted correlations 
 ###the normal scaled weights balances the weight adjustments for each pair of quantile.  this will likely help preserve normality in the data by re-scaling the weights in a balanced fashion. the products of two normal data sets will not necessarily be normal
  ##this method will iteratively update the scaled weights **only if*** they pass the shapiro test
  ###FIX ME: should consider scaling by the connected-ness of the modules.
  ###FIX ME: should consider iterating this processes across related pathways
  ###
   
 #  asymp.T<-sqrt(nrow(weightedAssociation)-2)*weightedAssociation/(2*sqrt(abs(1-weightedAssociation^2/var(ranking.weight)))) 
 asymp.T<-sqrt(ncol(weightedAssociation)-2)*weightedAssociation/sqrt(1-weightedAssociation^2)
 
  weighted.pvalue<- 2 * pt(abs(asymp.T), ncol(weightedAssociation) - 2, lower.tail = FALSE)

  xy<-data.frame(w=as.vector(weightedAssociation),t=as.vector(asymp.T),p=as.vector(weighted.pvalue),stringsAsFactors=FALSE)
 colR<-factor(c(rep("signf",length(which(xy$p<0.09))),rep("not.Signf", length(which(xy$p>=0.09)))))
  if(length(which(xy$p<0.09))>0){
 colR[which(xy$p<0.09)]<-factor("signf")
  }
 if(length(which(xy$p>=0.09))>0){
 colR[which(xy$p>=0.09)]<-factor("not.Signf")
  }
 xy2<-cbind(xy,colR)

 pp<-ggplot(data=xy2,aes(x=w,y=t,colour=colR))+geom_point()+geom_line(color='black',alpha=0.4)+ggtitle(paste0("Weighted Correlation Shapiro=",shapiro.test(weightedAssociation)$p.value," N=",length(names(df)) ) )
 
###plot normal fit
###the weighted association is a phase shift of a correlation values.
#weighted.fit<-fitdist(as.vector(weightedAssociation),"t",df=nSamples-2)
weighted.fit<-fitdist(as.vector(weightedAssociation),"norm")
plot(weighted.fit)
readkey()
n.sims<-20000
  stats <- replicate(n.sims, { r<-rnorm(n=length(weightedAssociation),
  mean=weighted.fit$estimate["mean"],
  sd=weighted.fit$estimate["sd"])
  as.numeric(ks.test(r,"pnorm",mean=weighted.fit$estimate["mean"],sd=weighted.fit$estimate["sd"])$statistic)})

  fit<-logspline(stats)
  weighted.sim.pvalue<-1-plogspline(ks.test(as.vector(weightedAssociation),
          "pnorm",
          mean=weighted.fit$estimate["mean"],
          sd=weighted.fit$estimate["sd"])$statistic,
          fit)
par(mfrow=c(1,1))
plot(ecdf(stats), las = 1, main = paste0("KS-test simulation Weighted p.value:",weighted.sim.pvalue), col = "darkorange", lwd = 1.7)
 readkey()
#####


#####normal unweighted call
  asymp.Cor.T<-sqrt(ncol(corrMap)-2)*corrMap/sqrt(1-corrMap^2)
  cor.pvalue<-2*pt(abs(asymp.Cor.T),ncol(corrMap)-2,lower.tail=FALSE)
 XY<-data.frame(w=as.vector(corrMap),t=as.vector(asymp.Cor.T),p=as.vector(cor.pvalue))
 colR2<-factor(c(rep("signf",length(which(XY$p<0.09))),rep("not.Signf", length(which(XY$p>=0.09)))))

 if(length(which(XY$p<0.09))>0){
 colR2[which(XY$p<0.09)]<-factor("signf")
  }
 if(length(which(XY$p>=0.09))>0){
 colR2[which(XY$p>=0.09)]<-factor("not.Signf")
 }
 XY2<-cbind(XY,colR2)

 pp2<-ggplot(data=XY2,aes(x=w,y=t,colour=colR2))+geom_point()+geom_line(color='black',alpha=0.4)+ggtitle(paste0("Unweighted Correlation Shapiro:",signif(shapiro.test(corrMap)$p.value,2) ))
  #######
####plot unweighted fit
unweighted.fit<-fitdist(as.vector(corrMap),"norm")
plot(unweighted.fit)
readkey()

  unweight.stats <- replicate(n.sims, { r<-rnorm(n=length(corrMap),
  mean=unweighted.fit$estimate["mean"],
  sd=unweighted.fit$estimate["sd"])
  as.numeric(ks.test(r,"pnorm",mean=unweighted.fit$estimate["mean"],sd=unweighted.fit$estimate["sd"])$statistic)})

  unweight.fit<-logspline(unweight.stats)
  unweighted.sim.pvalue<-1-plogspline(ks.test(as.vector(corrMap),
          "pnorm",
          mean=unweighted.fit$estimate["mean"],
          sd=unweighted.fit$estimate["sd"])$statistic,
          unweight.fit)
par(mfrow=c(1,1))
plot(ecdf(unweight.stats), las = 1, main = paste0("KS-test simulation Unweighted p.value:",unweighted.sim.pvalue), col = "darkorange", lwd = 1.7)
 readkey()



########
 ##winsorized
   win.weights<-ranking.weight
    x<-corrMap*win.weights
      large<-which(x>1)
      small<-which(x< ( -1) )
     wind<-x
     wind[large]<-1
     wind[small]<-(-1)
     

  asymp.W<-sqrt(ncol(wind)-2)*wind/(sqrt(1-wind^2))
   wind.pvalue<- 2 * pt(abs(asymp.W), ncol(wind) - 2, lower.tail = FALSE)
   wind.pvalue[is.na(wind.pvalue)]<-1
  wxy<-data.frame(w=as.vector(wind),t=as.vector(asymp.W),p=as.vector(wind.pvalue))
 colW<-factor(c(rep("signf",length(which(wxy$p<0.09))),rep("not.Signf", length(which(wxy$p>=0.09)))))

 if(length(which(wxy$p<0.09))>0){
 colW[which(wxy$p<0.09)]<-factor("signf")
  }
 if(length(which(wxy$p>=0.09))>0){
 colW[which(wxy$p>=0.09)]<-factor("not.Signf")
  }
  wxy2<-cbind(wxy,colW)

 ppw<-ggplot(data=wxy2,aes(x=w,y=t,colour=colW))+geom_point()+geom_line(color='black',alpha=0.4)+ggtitle(paste0("Windsorized Correlation Shapiro:",signif(shapiro.test(wind)$p.value,2) ))
######windsorized fit
windy.fit<-fitdist(as.vector(wind),"norm")
plot(windy.fit)
readkey()

 wind.stats <- replicate(n.sims, { r<-rnorm(n=length(wind),
  mean=windy.fit$estimate["mean"],
  sd=windy.fit$estimate["sd"])
  as.numeric(ks.test(r,"pnorm",mean=windy.fit$estimate["mean"],sd=windy.fit$estimate["sd"])$statistic)})

  win.fit<-logspline(wind.stats)
  wind.sim.pvalue<-1-plogspline(ks.test(as.vector(wind),
          "pnorm",
          mean=windy.fit$estimate["mean"],
          sd=windy.fit$estimate["sd"])$statistic,
          win.fit)
par(mfrow=c(1,1))
plot(ecdf(wind.stats), las = 1, main = paste0("KS-test simulation Windsorized p.value:",wind.sim.pvalue), col = "darkorange", lwd = 1.7)
 readkey()




############


  grid.arrange(pp,pp2,ppw,nrow=2,ncol=2)
  readkey()
   
 par(mfrow=c(2,2))
   hist(weighted.pvalue,main=paste0("Weighted Values Shapiro:",signif(shapiro.test(weightedAssociation)$p.value,2)))
   hist(corrMap.pvalue,main=paste0("Standard Unweighted Shapiro:",signif(shapiro.test(corrMap)$p.value,2)) )
   hist(wind.pvalue,main=paste0("Winsorized Shapiro:",signif(shapiro.test(wind)$p.value,2) ))
  readkey()


   pdf(paste0("Method_",how,"_Comparisons_Weighted_Correlations.pdf"))
  hist(weighted.pvalue,main=paste0("Weighted Values Shapiro:",signif(shapiro.test(weightedAssociation)$p.value,2)))
   hist(corrMap.pvalue,main=paste0("Standard Unweighted Shapiro:",signif(shapiro.test(corrMap)$p.value,2)) )
   hist(wind.pvalue,main=paste0("Winsorized Shapiro:",signif(shapiro.test(wind)$p.value,2) ))
  grid.arrange(pp,pp2,ppw,nrow=2,ncol=2)

   plot(weighted.fit)
   plot(unweighted.fit)
   plot(windy.fit)
  plot(ecdf(stats), las = 1, main = paste0("KS-test simulation Weighted p.value:",signif(weighted.sim.pvalue,1)), col = "darkorange", lwd = 1.7)
  plot(ecdf(unweight.stats), las = 1, main = paste0("KS-test simulation Unweighted p.value:",signif(unweighted.sim.pvalue,1)), col = "darkorange", lwd = 1.7)
  plot(ecdf(wind.stats), las = 1, main = paste0("KS-test simulation Windsorized p.value:",signif(wind.sim.pvalue,1)), col = "darkorange", lwd = 1.7)
 dev.off()
 ####################

  ha_mix_top=HeatmapAnnotation(density_line = anno_density(weightedAssociation,
                                            type = "line"),
                              heatmap = anno_density(weightedAssociation,
                                            type = "heatmap"),width=unit(4,"cm"))
  weighted.activation.direction<-Heatmap(asinh(activation.direction*ranking.weight),
 name="Activation Direction",
 col=colorRamp2(c(min(asinh(activation.direction)),0,max(asinh(activation.direction))),c("black","purple","yellow")),
 show_row_names=FALSE)
  
    weight<-Heatmap(asinh(weightedAssociation),
                   col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  name = "asinh(weight)",
                  top_annotation = ha_mix_top,
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=10),
                  column_title=paste0("Weighted Module ",how),
                  cell_fun=function(j,i,x,y,w,h,col){
                  asymp.T<-sqrt(ncol(weightedAssociation)-2)*weightedAssociation/(sqrt((1-weightedAssociation^2)))
                 weighted.pvalue<- 2 * pt(abs(asymp.T), ncol(weightedAssociation) - 2, lower.tail = FALSE)
                  if(weighted.pvalue[i,j]<0.09){
                   grid.text(sprintf("%.3f", weighted.pvalue[i,j]),x,y)
                  }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
                  })
  
      ###uses windsor method
      
   winsorized.heat<-Heatmap(asinh(wind),
                   col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  name = "asinh(weight)",
                  top_annotation = ha_mix_top,
                  top_annotation_height = unit(3, "cm"),
                  column_names_gp=gpar(fontsize=10),
                  column_title=paste0("Windsorized Weighted Module ",how),
                  cell_fun=function(j,i,x,y,w,h,col){
                  weighted.Pvalue<-corPvalueStudent(wind,ncol(wind)-2)
                   if(weighted.Pvalue[i,j]<0.06){
                   grid.text(sprintf("%.3f", weighted.Pvalue[i,j]),x,y)
                  }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
                  })
         
  
  draw(weight+weighted.activation.direction)
 readkey()
draw(winsorized.heat+weighted.activation.direction)
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
               asymp.Cor.T<-sqrt(ncol(corrMap)-2)*corrMap/sqrt(1-corrMap^2)
          cor.pvalue<-2*pt(abs(asymp.Cor.T),ncol(corrMap)-2,lower.tail=FALSE)
                 cor.pvalue[which(is.na(cor.pvalue))]<-1
                  if(cor.pvalue[i,j]<0.09){
                   grid.text(sprintf("%.3f", cor.pvalue[i,j]),x,y)
                  }
                   grid.rect(x,y,w,h,gp=gpar(fill=NA,col="black"))
                  })
   unweighted.activation.direction<-Heatmap(asinh(activation.direction),
        name="Activation Direction",
  col=colorRamp2(c(min(asinh(activation.direction)),0,max(asinh(activation.direction))),c("black","purple","yellow")),     
            show_row_names=FALSE)



    draw(unweight+unweighted.activation.direction)
   readkey()
 ###SO FAR SO GOOD>>> Fix rest.

 
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "");

 stopifnot(dim(textMatrix)==dim(moduleTraitCor))

 par(mfrow=c(1,1),mar = c(6, 10, 3, 3));
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
   draw(winsorized.heat+weighted.activation.direction)
    draw(unweight+unweighted.activation.direction)
   dev.off()
 
}
