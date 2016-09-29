#' @title this will query specific genes in the filteredCorSet lists
#' @description in wgcna_analysis a filteredCorSet list is constructed for trait gene correlations and module gene correlations, this will query for specific genes
wgcna_queryFilteredCorSet<-function(filteredCorSet,NKLigs=c("CLEC2B","KIR3DL3","CLEC4F","CLEC2A","KLRG2","CLEC4G","ULBP3","RAET1E","RAET1L","CLEC4C","CLEC2L","KLRD1","CLEC10A","CLEC3A","CLEC3B","CLEC19A","CLEC1A","MICA","RAET1G"),ISG=c("IFNGR1","TLR2","IFNAR2","IFITM3","IFNL4","IFNLR1","TLR9","IFNL1","CCL16","TRIL","CCL4","TLR5","IFNL2","IFNL3","CCL17","CCL22","TLR7","TLR10","IFITM5","IFNB1","CCL25","CCL15","IFRD2","ISG20"),ISG.pre=c("IL","DDX","IFN","TLR","IFIT","CCL","ISG","IRF","NFKB"),NK.pre=c("CLEC","KIR","KLR","RAET","MIC","CLEX","ULBP")){

trait.qFilter<-filteredCorSet[["trait.qFilter"]]
module.qFilter<-filteredCorSet[["module.qFilter"]]

  for(i in 1:length(trait.qFilter)){
    traitidx<-vector()
    for(j in 1:length(ISG.pre)){
     id<-grep(ISG.pre[j],trait.qFilter[[i]]$hgnc_symbol) 
      if(length(id)>0){
         traitidx<-c(traitidx,id)  
         isg<-trait.qFilter[[i]][traitidx,]
         write.csv(isg,file=paste0(names(trait.qFilter)[i],"-ISGQuery.csv"),quote=FALSE,row.names=TRUE)
        save(isg,file=paste0(names(trait.qFilter)[i],"-ISGQuery.RData"),compress=TRUE)
      }
   } #for j
  } #for i 

for(i in 1:length(trait.qFilter)){
    traitidx<-vector()
    for(j in 1:length(NK.pre)){
     id<-grep(NK.pre[j],trait.qFilter[[i]]$hgnc_symbol)
      if(length(id)>0){
         traitidx<-c(traitidx,id)
         nk<-trait.qFilter[[i]][traitidx,]
         write.csv(nk,file=paste0(names(trait.qFilter)[i],"-NKQuery.csv"),quote=FALSE,row.names=TRUE)
        save(nk,file=paste0(names(trait.qFilter)[i],"-NKQuery.RData"),compress=TRUE)
      }
   } #for j
  } #for i 
#########################################
###FIX ME: find which modules are associated with ISG or NK

### each module gene expression is up or down ?
### each module gene enrichment positive or negative? this gives direction to the direction of neg. or pos. correlation
########################################

 for(i in 1:length(module.qFilter)){
  traitidx<-vector()
   df<-module.qFilter[[i]]
   id<-which(abs(df[1,])>=0.5)
   df[id,] 
   readkey()
}

} #main
