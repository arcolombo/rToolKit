#' @title filters low mean counts 
#' @description often it is needed to filter counts, this does that
#' @param datExpr from the wgcna object lnames
#' @param annot from the gene annotation call
#' @export
#' @import WGCNA
#' @return a data frame with p adjusted q value genes
wgcna_filter<-function(geneTraitModuleDF,datExpr,annot){


##the geneTraitModuleDF has pValues, there are not Fisher pvalues....  maybe add ?
##has a qvalue for each module color
q.list<-lapply(geneTraitModuleDF[,grepl("p.MM",colnames(geneTraitModuleDF))],
                function(x) qvalue(x)$qvalues)
names(q.list)<-sub("^p.","q.",names(q.list))

###module list p and q value
colors<-sub("^MM","",colnames(geneTraitModuleDF)[grepl("^MM",colnames(geneTraitModuleDF))])
###need to subset each individual module color genes,corrl,p.value,fisherp,qvalue  in a list.  filter and then write out to text files each module

#######################################


meanExpressionByArray=apply( datExpr,1,mean, na.rm=T)  
NumberMissingByArray=apply( is.na(data.frame(datExpr)),1, sum)  
all(NumberMissingByArray==0)

  y1<-length(grep("pHSC_",rownames(datExpr)))
  y2<-length(grep("LSC_",rownames(datExpr)))
  y3<-length(grep("Blast_",rownames(datExpr)))
  y<-c(rep(1,y1),rep(2,y2),rep(3,y3))

  GS1= as.numeric(cor(y, datExpr, use="p"))
  # Network terminology: GS1 will be referred to as signed gene significance measure
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  p.Standard=corPvalueFisher(GS1, nSamples =length(y) )
  # since the q-value function has problems with missing data, we use the following trick
  p.Standard2=p.Standard
  p.Standard2[is.na(p.Standard)]=1
  q.Standard=qvalue(p.Standard2)$qvalues
  # Form a data frame to hold the results
  GeneName=colnames(datExpr)
  hugoID<-match(GeneName,annot$ensembl_gene_id)
  HugoName<-annot$hgnc_symbol[hugoID]
  StandardGeneScreeningResults=data.frame(GeneName,
                               PearsonCorrelation=GS1, 
                               p.Standard,
                               q.Standard,
                               HGNC=HugoName)
message("FDR under q=0.20")
print(table(q.Standard<.20))
return(StandardGeneScreeningResults)
} 
