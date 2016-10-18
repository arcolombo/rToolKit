#' @title this compares genes from across two different experiments
#' @description comparing genes across different platforms to compare and contrast gene networks across experiments ran on different platforms.
compareGenesAcross<-function(kexp1,kexp2,whichWGCNA="single",species="Homo.sapiens",read.cutoff=2,softPower=10,how="cpm"){
  #FIX ME: Add TPM Support
  ##FIX ME: break parts of this into sub=functions
  ##FIX ME:  finish the Miller tutorial and then build the entire data as a full data base and re-run the same analysis on the full data set, not just the first block.  
   ##FIX ME: Need Network Statistics. Horvath Text
  ##TO DO:  need to run blockwise on ens80 buenrostro and blockwise on tcga
  ##then this will just find the common genes and begin the analysis given common genes and common repeats found between the two.  
  ##either 1: assign modules 1 to tcga, or merely overlap them 
  ###I do not fully unerstand the ME and kMEs yet.  
 
  whichWGCNA<-match.arg(whichWGCNA,c("single","block"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  how<-match.arg(how,c("cpm","tpm"))
  if(how=="cpm"){
  ###rows must be SAMPLES columns genes
  cpm1<-collapseBundles(kexp1,"gene_id",read.cutoff=read.cutoff)
  cpm1<-cpm1[!grepl("^ERCC",rownames(cpm1)),]
  cpm1<-cpm1[grepl("^ENSG",rownames(cpm1)),]
 ###rows must be SAMPLES columns genes
  cpm2<-collapseBundles(kexp2,"gene_id",read.cutoff=read.cutoff)
  cpm2<-cpm2[!grepl("^ERCC",rownames(cpm2)),]
  cpm2<-cpm2[grepl("^ENSG",rownames(cpm2)),]
  d1<-DGEList(counts=cpm1)
  cpm.norm1<-cpm(d1,normalized.lib.sizes=TRUE)
  cpm1<-cpm.norm1
  cpm.norm1<-NULL
  d2<-DGEList(counts=cpm2)
  cpm.norm2<-cpm(d2,normalized.lib.sizes=TRUE)
  cpm2<-cpm.norm2
  cpm.norm2<-NULL
 match.id<-rownames(cpm2)%in%rownames(cpm1)
 cpm22<-cpm2[match.id,]
 m.id<-match(rownames(cpm22),rownames(cpm1))  
 cpm1<-cpm1[m.id,]
 stopifnot(all(rownames(cpm22)==rownames(cpm1))==TRUE)
 cpm2<-cpm22

 ##grabbing repeats
  ###rows must be SAMPLES columns genes
  rexp1<-findRepeats(kexp1)
  rpm1<-collapseBundles(rexp1,"gene_id",read.cutoff=read.cutoff)
  rpm1<-rpm1[!grepl("^ERCC",rownames(rpm1)),]
  ###rows must be SAMPLES columns genes
  rexp2<-findRepeats(kexp2)
  rpm2<-collapseBundles(rexp2,"gene_id",read.cutoff=read.cutoff)
  rpm2<-rpm2[!grepl("^ERCC",rownames(rpm2)),]
  
  rd1<-DGEList(counts=rpm1)
  rpm.norm1<-cpm(rd1,normalized.lib.sizes=TRUE)
  rpm1<-rpm.norm1
  rpm.norm1<-NULL
  rd2<-DGEList(counts=rpm2)
  rpm.norm2<-cpm(rd2,normalized.lib.sizes=TRUE)
  rpm2<-rpm.norm2
  rpm.norm2<-NULL
 rmatch.id<-rownames(rpm2)%in%rownames(rpm1)
 rpm22<-rpm2[rmatch.id,]
 rm.id<-match(rownames(rpm22),rownames(rpm1))
 rpm1<-rpm1[rm.id,]
 stopifnot(all(rownames(rpm22)==rownames(rpm1))==TRUE)
 rpm2<-rpm22

 } else{
 stop("Not supporting TPMs yet...")
 }
 

 datExprG1<-cpm1 ##aml genes
 datExprG2<-cpm2 ##TCGA genes
 datExprR1<-rpm1 ##AML repeats
 datExprR2<-rpm2 ##TCGA repeats
 ##rank the row means
  rankExprG1<-rank(rowMeans(datExprG1))
  rankExprG2<-rank(rowMeans(datExprG2))
  random5000<-sample(rownames(datExprG1),5000)
  rankConnG1<-rank(softConnectivity(t(datExprG1[random5000,] ),type="signed",power=softPower))
  rankConnG2<-rank(softConnectivity(t(datExprG2[random5000,] ),type="signed",power=softPower))

 
  rankExprR1<-rank(rowMeans(datExprR1))
  rankExprR2<-rank(rowMeans(datExprR2))
  random400<-sample(rownames(datExprR1),400)
  rankConnR1<-rank(softConnectivity(t(datExprR1[random400,]),type="signed",power=softPower))
  rankConnR2<-rank(softConnectivity(t(datExprR2[random400,]),type="signed",power=softPower))


 #here we correlate ranked average gene expression and overall connectivity. high correlation means we can find similarities between sets.
 par(mfrow=c(2,2))
  verboseScatterplot(rankExprG1,rankExprG2, xlab="Ranked Expression (G1)",
  ylab="Ranked Expression (G2)")
  verboseScatterplot(rankConnG1,rankConnG2, xlab="Ranked Connectivity (G1)",
  ylab="Ranked Connectivity (G2)")
  verboseScatterplot(rankExprR1,rankExprR2, xlab="Ranked Expression (R1)",
  ylab="Ranked Expression (R2)")
  verboseScatterplot(rankConnR1,rankConnR2, xlab="Ranked Connectivity (R1)",
  ylab="Ranked Connectivity (R2)") 
 readkey()
 pdf("GeneralNetworkProperies.pdf",width=9,height=10)
par(mfrow=c(2,2))
  verboseScatterplot(rankExprG1,rankExprG2, xlab="Ranked Expression (G1)",
  ylab="Ranked Expression (G2)")
  verboseScatterplot(rankConnG1,rankConnG2, xlab="Ranked Connectivity (G1)",
  ylab="Ranked Connectivity (G2)")
  verboseScatterplot(rankExprR1,rankExprR2, xlab="Ranked Expression (R1)",
  ylab="Ranked Expression (R2)")
  verboseScatterplot(rankConnR1,rankConnR2, xlab="Ranked Connectivity (R1)",
  ylab="Ranked Connectivity (R2)")
 dev.off()
 
 ##we need a cluster computer to compute all genes, so we take the first 5000 genes with the highest expression.  to get the full set HPCC is required.
  keepGenesExpr = rank(-rowMeans(datExprG1))<=5000
datExprg1 = datExprG1[keepGenesExpr,]
datExprg2 = datExprG2[keepGenesExpr,]
  ###FIX ME:::: add this for all genes blockwise.
 ##create adjacency matrices
  adjacencyG1 = adjacency(t(datExprg1),power=softPower,type="signed");
  diag(adjacencyG1)=0
  dissTOMG1 = 1-TOMsimilarity(adjacencyG1, TOMType="signed")
  geneTreeG1 = flashClust(as.dist(dissTOMG1), method="average")
  adjacencyG2 = adjacency(t(datExprg2),power=softPower,type="signed");
  diag(adjacencyG2)=0
  dissTOMG2 = 1-TOMsimilarity(adjacencyG2, TOMType="signed")
  geneTreeG2 = flashClust(as.dist(dissTOMG2), method="average") 
  ##for Repeats
  adjacencyR1 = adjacency(t(datExprR1),power=softPower,type="signed");
  diag(adjacencyR1)=0
  dissTOMR1 = 1-TOMsimilarity(adjacencyR1, TOMType="signed")
  geneTreeR1 = flashClust(as.dist(dissTOMR1), method="average")
  adjacencyR2 = adjacency(t(datExprR2),power=softPower,type="signed");
  diag(adjacencyR2)=0
  dissTOMR2 = 1-TOMsimilarity(adjacencyR2, TOMType="signed")
  geneTreeR2 = flashClust(as.dist(dissTOMR2), method="average") 
 ##now to compare networks
  par(mfrow=c(2,2))
  plot(geneTreeG1,
       xlab="",
       sub="",
       main="Gene clustering on TOM-based dissimilarity (G1)",
       labels=FALSE,hang=0.04);
  plot(geneTreeG2,
       xlab="",
       sub="",
       main="Gene clustering on TOM-based dissimilarity (G2)",
       labels=FALSE,hang=0.04);
  plot(geneTreeR1,
       xlab="",sub="",
       main="Gene clustering on TOM-based dissimilarity (R1)",
      labels=FALSE,hang=0.04);
  plot(geneTreeR2,
       xlab="",
       sub="",
       main="Gene clustering on TOM-based dissimilarity (R2)",
       labels=FALSE,hang=0.04);
  readkey()
 pdf("dendrogram.pdf",height=8,width=16)
 par(mfrow=c(2,2))
  plot(geneTreeG1,
     xlab="",
     sub="",
     main="Gene clustering on TOM-based dissimilarity (A1)",
     labels=FALSE,hang=0.04);
  plot(geneTreeG2,
     xlab="",
     sub="",
     main="Gene clustering on TOM-based dissimilarity (A2)",
     labels=FALSE,hang=0.04);
 plot(geneTreeR1,
       xlab="",sub="",
       main="Gene clustering on TOM-based dissimilarity (R1)",
      labels=FALSE,hang=0.04);
  plot(geneTreeR2,
       xlab="",
       sub="",
       main="Gene clustering on TOM-based dissimilarity (R2)",
       labels=FALSE,hang=0.04);
  dev.off()

  ##it is desirable to select the module size, in some case it is better to have fewer modules with more genes.and vice versa.
  mColorG1=NULL
  mColorG2=NULL
  mColorR1=NULL
  mColorR2=NULL
  for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeG1, pamStage=FALSE,
 minClusterSize = (30-3*ds), cutHeight = 0.99,
 deepSplit = ds, distM = dissTOMG1)
 mColorG1=cbind(mColorG1,labels2colors(tree$labels));
  }
 for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeG2, pamStage=FALSE,
 minClusterSize = (30-3*ds), cutHeight = 0.99,
 deepSplit = ds, distM = dissTOMG2)
 mColorG2=cbind(mColorG2,labels2colors(tree$labels));
  }
 for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeR1, pamStage=FALSE,
 minClusterSize = (30-3*ds), cutHeight = 0.99,
 deepSplit = ds, distM = dissTOMR1)
 mColorR1=cbind(mColorR1,labels2colors(tree$labels));
  }
 for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeR2, pamStage=FALSE,
 minClusterSize = (30-3*ds), cutHeight = 0.99,
 deepSplit = ds, distM = dissTOMR2)
 mColorR2=cbind(mColorR2,labels2colors(tree$labels));
  }
  plotDendroAndColors(geneTreeG1, mColorG1, paste("dpSplt =",0:3),dendroLabels=FALSE,main="Top G1 Modules Split Level");
  deepSplit<-readDeepSplit()
  modulesG1 = mColorG1[,deepSplit] 
  readkey()
  ##
   plotDendroAndColors(geneTreeG2, mColorG2, paste("dpSplt =",0:3), main = "Top G2 Modules Split Level",dendroLabels=FALSE);
  deepSplit<-readDeepSplit()
  modulesG2 = mColorG2[,deepSplit]
  readkey()
  ##
 plotDendroAndColors(geneTreeR1, mColorR1, paste("dpSplt =",0:3), main = "Repeats R1 Module Split Level",dendroLabels=FALSE);
  deepSplit<-readDeepSplit()
  modulesR1 = mColorR1[,deepSplit]
  readkey()
  ####
 plotDendroAndColors(geneTreeR2, mColorR2, paste("dpSplt =",0:3), main = "Repeats R2 Module Split Level",dendroLabels=FALSE);
  deepSplit<-readDeepSplit()
  modulesR2 = mColorR2[,deepSplit]
  readkey()

  pdf("Module_choices.pdf", height=10,width=25);
   par(mfrow=c(2,2))
   plotDendroAndColors(geneTreeG1, mColorG1, paste("dpSplt =",0:3),dendroLabels=FALSE,main="Top G1 Modules Split Level");
   ##
   plotDendroAndColors(geneTreeG2, mColorG2, paste("dpSplt =",0:3), main = "Top G2 Modules Split Level",dendroLabels=FALSE);
  ##
 plotDendroAndColors(geneTreeR1, mColorR1, paste("dpSplt =",0:3), main = "Repeats R1 Module Split Level",dendroLabels=FALSE);
  ####
 plotDendroAndColors(geneTreeR2, mColorR2, paste("dpSplt =",0:3), main = "Repeats R2 Module Split Level",dendroLabels=FALSE);
   dev.off()


  ##calculate PCs and bar graphs
  PCsG1 = moduleEigengenes(t(datExprg1), colors=modulesG1)
  ME_G1 = PCsG1$eigengenes
  distPCG1 = 1-abs(cor(ME_G1,use="p"))
  distPCG1 = ifelse(is.na(distPCG1), 0, distPCG1)
  pcTreeG1 = hclust(as.dist(distPCG1),method="a")
  MDS_G1 = cmdscale(as.dist(distPCG1),2)
  colorsG1 = names(table(modulesG1))
 

  PCsG2 = moduleEigengenes(t(datExprg2), colors=modulesG2)
  ME_G2 = PCsG2$eigengenes
  distPCG2 = 1-abs(cor(ME_G2,use="p"))
  distPCG2 = ifelse(is.na(distPCG2), 0, distPCG2)
  pcTreeG2 = hclust(as.dist(distPCG2),method="a")
  MDS_G2 = cmdscale(as.dist(distPCG2),2)
  colorsG2 = names(table(modulesG2))

  PCsR1 = moduleEigengenes(t(datExprR1), colors=modulesR1)
  ME_R1 = PCsR1$eigengenes
  distPCR1 = 1-abs(cor(ME_R1,use="p"))
  distPCR1 = ifelse(is.na(distPCR1), 0, distPCR1)
  pcTreeR1 = hclust(as.dist(distPCR1),method="a")
  MDS_R1 = cmdscale(as.dist(distPCR1),2)
  colorsR1 = names(table(modulesR1))

  PCsR2 = moduleEigengenes(t(datExprR2), colors=modulesR2)
  ME_R2 = PCsR2$eigengenes
  distPCR2 = 1-abs(cor(ME_R2,use="p"))
  distPCR2 = ifelse(is.na(distPCR2), 0, distPCR2)
  pcTreeR2 = hclust(as.dist(distPCR2),method="a")
  MDS_R2 = cmdscale(as.dist(distPCR2),2)
  colorsR2 = names(table(modulesR2))

 
       
        
 moduleVisualization(pcTree=pcTreeG1,
                     MDS_G=MDS_G1,
                     colorsG=colorsG1,
                     treeOrder=geneTreeG1$order,
                     modulesG=modulesG1,
                     moduleEigen=ME_G1,
                     imageName="G1")

       
 moduleVisualization(pcTree=pcTreeG2,
                     MDS_G=MDS_G2,
                     colorsG=colorsG2,
                     treeOrder=geneTreeG2$order,
                     modulesG=modulesG2,
                     moduleEigen=ME_G2,
                     imageName="G2")

        
 moduleVisualization(pcTree=pcTreeR1,
                     MDS_G=MDS_R1,
                     colorsG=colorsR1,
                     treeOrder=geneTreeR1$order,
                     modulesG=modulesR1,
                     moduleEigen=ME_R1,
                     imageName="R1")


moduleVisualization(pcTree=pcTreeR2,
                     MDS_G=MDS_R2,
                     colorsG=colorsR2,
                     treeOrder=geneTreeR2$order,
                     modulesG=modulesR2,
                     moduleEigen=ME_R2,
                     imageName="R2")

 ####Now to compare across kexp1 and kexp2
 preserve<-analyzePreservation(geneTree1=geneTreeG1,
                               module1=modulesG1,
                               main1="Gene Dendrogram G1 and Modules G1",
                               geneTree2=geneTreeG2,
                               main2="Gene Dendrogram G2 and modules G1")

  preservationScore<-quantifyPreservation(datExpr1=datExprg1,
                                          datExpr2=datExprg2,
                                          modules=modulesG1)
    readkey()
  repeat.preserve<-analyzePreservation(geneTree1=geneTreeR1,
                               module1=modulesR1,
                               main1="Gene Dendrogram R1 and Modules R1",
                               geneTree2=geneTreeR2,
                               main2="Gene Dendrogram R2 and modules R1")

  repeat.preservationScore<-quantifyPreservation(datExpr1=datExprR1,
                                          datExpr2=datExprR2,
                                          modules=modulesR1)
 
   ###calculate connectivity of ME.  kME measures correlatoin of gene and Module Eigenvalues. 
  #kME
  genekME<-analyzekME(datExpr1=datExprg1,
                      datExpr2=datExprg2, 
                      ME_G1=ME_G1,
                      colorsG1=colorsG1,
                      modulesG1=modulesG1,
                      main="Gene")
 repeatkME<-analyzekME(datExpr1=datExprR1,
                      datExpr2=datExprR2,
                      ME_G1=ME_R1,
                      colorsG1=colorsR1,
                      modulesG1=modulesR1,
                      main="Repeat")

} #Main
