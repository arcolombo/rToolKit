#' @title this analyzes the filtered correlations
#' @import WGCNA
#' @import cluster
wgcna_filteredAnalysis<-function(kexp,lnames,useBiCor=TRUE){

 bwModuleColors<-lnames[["moduleColors"]]
  MEs<-lnames[["MEs"]]
  datExpr<-lnames[["datExpr"]]
  datTraits<-lnames[["datTraits"]]
  annot<-lnames[["annot"]]
  MEs<-lnames[["MEs"]]
  moduleTraitCor<-lnames[["moduleTraitCor"]]
  moduleTraitPvalue<-lnames[["moduleTraitPvalue"]]
  modNames = substring(names(MEs), 3)
  nGenes<-ncol(datExpr)
  nSamples = nrow(datExpr);

  if(useBiCor==TRUE){
  ADJ1=abs(bicor(datExpr,use="all.obs"))^6
  } else{
  ADJ1=abs(cor(datExpr,use="p"))^6
   }
# When you have relatively few genes (<5000) use the following code
  if(nGenes<5000){
  k=as.vector(apply(ADJ1,2,sum, na.rm=T))
  } else if(nGenes>5000){
# When you have a lot of genes use the following code
  k=softConnectivity(datE=datExpr,power=6) 
  }
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
datExpr=datExpr[, rank(-k,ties.method="first" )<=3600]
dissADJ=1-ADJ1

dissTOM=TOMdist(ADJ1)
collectGarbage()
pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)



pamTOM4=pam(as.dist(dissTOM), 4)
pamTOM5=pam(as.dist(dissTOM), 5)
pamTOM6=pam(as.dist(dissTOM), 6)





hierADJ=hclust(as.dist(dissADJ), method="average" )

colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")




branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )



colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ, 
                              cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

# Plot results of all module detection methods together:
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ, 
                    colors=data.frame(colorStaticADJ, 
                                     colorDynamicADJ, colorDynamicHybridADJ), 
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")





# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1 
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                       deepSplit=2, pamRespectsDendro = FALSE))
# Now we plot the results
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM, 
               colors=data.frame( colorStaticTOM, 
                                 colorDynamicTOM, colorDynamicHybridTOM), 
               dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
               main = "Gene dendrogram and module colors, TOM dissimilarity")




tabStaticADJ=table(colorStaticADJ)
tabStaticTOM=table(colorStaticTOM)
tabDynamicADJ=table(colorDynamicADJ)
tabDynamicTOM=table(colorDynamicTOM)
tabDynamicHybridADJ =table(colorDynamicHybridADJ)
tabDynamicHybridTOM =table(colorDynamicHybridTOM)

print(tabStaticADJ)
print(tabStaticTOM)
print(tabDynamicADJ)
print(tabDynamicTOM)
print(tabDynamicHybridADJ)
print(tabDynamicHybridTOM)


colorh1= colorDynamicHybridTOM

rm(ADJ1); rm(dissADJ);              
collectGarbage()

if(useBiCor==FALSE){
datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")



} else {
datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(bicor(datME, use="all.obs"), 2)
dissimME=(1-t(bicor(datME, use="all.obs")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
}

plotMEpairs(datME)

###plotting network from colorh1 as hybrid dynamic tree
cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1),  main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")



###colorh1 as colorDynamicTOM
power=6
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "TOM heatmap plot, module genes" )  ###slow



#########fix me : what is NS1 ?

##y is a clinical trait with one value per gene
NS1=networkScreening(y=y, datME=datME, datExpr=datExpr,
         oddPower=3, blockSize=1000, minimumSampleSize=4,
         addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
 
sizeGrWindow(7,7)
topList=rank(NS1$p.Weighted,ties.method="first")<=30
gene.names= names(datExpr)[topList]
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=FALSE,
                   power=1, main="signed correlations")




sizeGrWindow(7,7)
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")





sizeGrWindow(7,7)
# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=TRUE,
                   power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="D. TOM in an unsigned network")



}# main
