#' @title runs qusage to test enrichment of ENSEMBL IDs
#' @description rus qusage to find enriched sets
#' @import knitr
#' @import qusage
#' @import arkas
#' @import TxDbLite
#' @export
qusageAnalysis<-function(annotatedKexp,byLevel=c("transcript","gene"),reactomeSpecies=c("Mus musculus","Homo sapiens")){

reactomeSpecies<-match.arg(reactomeSpecies,c("Mus musculus","Homo sapiens"))
design<-metadata(annotatedKexp)$design
comparisonSamples<-rownames(design)[which(design[,2]>0)]
controlSamples<-rownames(design)[which(design[,2]<1)]

#Run enrichment analysis
if(byLevel=="transcript"){
#transcript level

trnxMap<-mapToReactome(rownames(annotatedKexp),
                       type="transcript",
                       species=reactomeSpecies,
                       build=84)

rSets<-reactomeSets(species=reactomeSpecies,type="transcript",mappedReactome=trnxMap)
idx<-which(rowSums(counts(annotatedKexp))>0)

filteredTrnx<-counts(annotatedKexp)[idx,]

tL<-log2(filteredTrnx+0.001)


 for(i in 1:length(comparisonSamples) ) {
 idx<-which(colnames(tL)==comparisonSamples[i])
 colnames(tL)[idx]<-paste0("COMP_0")
 }

   for(i in 1:length(controlSamples)) {
  idx2<-which(colnames(tL)==controlSamples[i])
  colnames(tL)[idx2]<-paste0("CNTRL_1")
   }



sampleSplits<-split(tL,colnames(tL))
sampleVars<-lapply(sampleSplits,var)
print("Sample Group Variances")
print(sampleVars)

#running speedSage
tx.Results<-qusage(tL,
                      colnames(tL),
                      "COMP_0-CNTRL_1",
                       rSets,
                       var.equal=FALSE,
                       n.points=2^12)

#print reactome plots
kable(summary(tx.Results),caption="Transcript Enrichment")
p.vals<-pdf.pVal(tx.Results)
q.vals<-p.adjust(p.vals,method="fdr")
tx.Stats<-data.frame(names(tx.Results$pathways),q.vals)
write.table(tx.Stats,file=paste0("transcript-enrichment-QValues.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

qs<-qsTable(tx.Results,number=numPathways(tx.Results))
if(numPathways(tx.Results) <=10) {
message("plotting confidence intervals of pathways")
plotCIs(tx.Results)
}

if(numPathways(tx.Results)>10) {
message(paste0("plotting ",numPathways(tx.Results), " detected pathways"))
plot(tx.Results)
}

#plot URL
message("printing Reactome Urls for exploration...")
tx.Url<- getReactomeUrl(qs[,1])
tx.DF<-data.frame(qs,tx.Url, stringsAsFactors=FALSE)

data(reactomePathways,package="TxDbLite")

indX<-names(reactomePathways) %in% tx.DF$pathway.name
for(i in 1:nrow(tx.DF)) {
inner<-which(tx.DF$pathway.name[i] == names(reactomePathways[indX]))
tx.DF$Pathway.Description[i]<-reactomePathways[indX][[inner]]
}
tx.DF$Pathway.Description<-format(tx.DF$Pathway.Description,justify="left")


trnxReactomePathways<-reactomePathways[indX]
kable(tx.DF,caption="Full Enrichment Transcript Analysis")
#save output

#printing a list of reactomeIDs per gene-name
rID.trnx<-lapply(rSets,function(x) mapHugo(x,byType="transcript"))
rID.trnx<-rID.trnx[lapply(rID.trnx,length)>0]

print("printing out the reactome pathway transcript sets and associated gene names")
print(rID.trnx)

save(rID.trnx,file=paste0("rID.trnx.RData"))
save(trnxReactomePathways,file=paste0("trnxReactomePathways.RData"))
save(rSets,file=paste0("rSets.RData"))
write.table(tx.DF,file=paste0("transcript-enrichment.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
save(tx.Results,file=paste0("tx.Results.RData"))
} #tx level

if(byLevel=="gene"){
#gene level
coll<-collapseBundles(annotatedKexp,"gene_id")

gnMap<-mapToReactome(rownames(coll),
                       type="gene",
                       species=reactomeSpecies,
                       build=84)

grSets<-reactomeSets(species=reactomeSpecies,type="gene",mappedReactome=gnMap)
gidx<-which(rowSums(coll)>0)
filteredGn<-coll[gidx,]
gtL<-log2(filteredGn+0.001)

if(any(is.infinite(gtL))==TRUE) {
message("I detected infinite values...")
}



 for(i in 1:length(comparisonSamples) ) {
 idx<-which(colnames(gtL)==comparisonSamples[i])
 colnames(gtL)[idx]<-paste0("COMP_0")
 }

   for(i in 1:length(controlSamples)) {
  idx2<-which(colnames(gtL)==controlSamples[i])
  colnames(gtL)[idx2]<-paste0("CNTRL_1")
   }



sampleSplits<-split(gtL,colnames(gtL))
sampleVars<-lapply(sampleSplits,var)
print("Sample Group Variances")
print(sampleVars)

gn.Results<-qusage(gtL,
                      colnames(gtL),
                      "COMP_0-CNTRL_1",
                       grSets,
                       var.equal=FALSE,
                       n.points=2^12)



#plot results
kable(summary(gn.Results),caption="Gene Enrichment Statistics")
p.vals<-pdf.pVal(gn.Results)
q.vals<-p.adjust(p.vals,method="fdr")

gn.Stats<-data.frame(names(gn.Results$pathways),q.vals)
write.table(gn.Stats,file=paste0("gene-enrichment-QValues.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

gs<-qsTable(gn.Results,number=numPathways(gn.Results))

if(numPathways(gn.Results) <=10) {
message("plotting confidence intervals of pathways")
plotCIs(gn.Results)
}

if(numPathways(gn.Results)>10) {
message(paste0("plotting ",numPathways(gn.Results), " detected pathways"))
plot(gn.Results)
}

#plot URLs

gn.Url<- getReactomeUrl(gs[,1])
gn.DF<-data.frame(gs,gn.Url, stringsAsFactors=FALSE)


data(reactomePathways,package="TxDbLite")

gindX<-names(reactomePathways) %in% gn.DF$pathway.name
for(i in 1:nrow(gn.DF)) {
ginner<-which(gn.DF$pathway.name[i] == names(reactomePathways[gindX]))
gn.DF$Pathway.Description[i]<-reactomePathways[gindX][[ginner]]
}
 gn.DF$Pathway.Description<-format(gn.DF$Pathway.Description,justify="left")

kable(gn.DF,caption="Full Enrichment Gene Analysis")

#printing a list of reactome IDs per hugo gene name
reactomeIDgeneID<-lapply(grSets,function(x) mapHugo(x,byType="gene"))


id<-gn.DF$pathway.name[(gn.DF$p.Value<0.1)]
print("the following geneSets with a p.value <0.1 have the following gene names in the set:")
print(reactomeIDgeneID[id])

save(reactomeIDgeneID,file=paste0("reactomeIDgeneID.RData"))

gnReactomePathways<-reactomePathways[gindX]


#save output
save(gn.Results,file=paste0("gn.Results.RData"))
save(gnReactomePathways,file=paste0("gnReactomePathways.RData"))
save(grSets,file=paste0("grSets.RData"))
write.table(gn.DF,file=paste0("gene-enrichment.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}

}#main
