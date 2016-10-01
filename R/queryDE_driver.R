#' @title queries a Differential Expression gene list for drivers in a given module biotype association
#' @description given a Differential Expression gene list it is a great interest to query the DE list for each driver identity.
#' @export
#' @return data print out to screen
queryDE_driver<-function(DE.list=NULL, moduleColor="blue", module.biotype.driver=c("Alu","Endogenous Retrovirus","ERV3","ERVK","ERVL","L1","L2"),path="~/Documents/Arkas-Paper-Data/AML-bonemarrow-LSCs/patient-plot-data/wgcna_data/repeat_short_list_biotypes/enrichment/RObjects/",p.cutoff=0.05){
  stopifnot(is.null(DE.list)==FALSE) 

###do this as a loop across biotypes for a fixed color
for(i in 1:length(module.biotype.driver)){

 if( file.exists(paste0(path,"/",moduleColor,"_GCor.",module.biotype.driver[i],".RData"))==TRUE){
   cat(paste0("found ",paste0(path,"/",moduleColor,"_GCor.",module.biotype.driver[i],".RData")))
  load(paste0(path,"/",moduleColor,"_GCor.",module.biotype.driver[i],".RData"))
   print(head(df))
  } else{
   stop(paste0("Did not find ",moduleColor," ",module.biotype.driver[i]))
   }
  potential<-df[which(df$p.value<p.cutoff),]
  cat(paste0("there are ",nrow(potential)," significant correlated repeats \n"))
  drivers<-DE.list[rownames(DE.list)%in%rownames(potential),]
  up<-drivers[which(drivers$logFC>0),]
  cat(paste0("there are ",nrow(up)," +logFC DE \n"))
  down<-drivers[which(drivers$logFC<0),]
  cat(paste0("thre are ",nrow(down)," -logFC DE \n"))
  stopifnot((nrow(up)+nrow(down))==nrow(drivers))
  
 ###test if down logFC is down correlated
  potential.down<-potential[which(potential[,2]<0),]
 cat(paste0("high significant -cor :",nrow(potential.down),"\n" ))
  potential.up<-potential[which(potential[,2]>0),]
 cat(paste0("high signif. +cor :",nrow(potential.up),"\n"))
  cat(paste0("for p.cutoff ",p.cutoff," most significant R^2XR^2 in DE.list: ",moduleColor,"X",module.biotype.driver[i]," ",100*(nrow(drivers)/nrow(potential)),"% constitute ",100*(nrow(drivers)/nrow(DE.list)),"%", "\n"))
  
 cat(paste0("all +logFC in DE correspond to +cor in ",moduleColor," ",module.biotype.driver[i],": ",all(rownames(up)%in%rownames(potential.up)),"\n" ))
 cat(paste0("all -logFC in DE correspond to -cor in ",moduleColor," ",module.biotype.driver[i],": ", all(rownames(down)%in%rownames(potential.down)),"\n" ) )
 # cat(paste0("all -logFC in DE correspond to +cor in ",moduleColor," ",module.biotype.driver[i],": ", all(rownames(down)%in%rownames(potential.down)),"\n" ) )
  

readkey()
 
} #for loop mods
###fix me: add a print out report


##the potential drivers that are up, queryDF  and down queryDF


##FIX ME: Must add a consensus, perform a set analysis on each potential.down potential.up.  what potential.up/down are common among the biotypes for a module color?

 return(NULL)
} #main_
