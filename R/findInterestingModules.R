#' read key for inputing interesting modules
#' @export
findInterestingModules<-function(n.modules=7){

   intModules<-matrix(nrow=n.modules,ncol=1)
   message(paste0("please select ",n.modules," interesting modules"))
   for(i in 1:n.modules){
   intModules[i,1]<-.readModules(i)
    }
   return(as.vector(intModules))
} #main


.readModules<-function(i){
   cat("Enter [color] only :\n")
    line <- readline()
    message(paste0("Color Selected ",i," :",line))
    return(as.character(line))
}


