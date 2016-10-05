#' @title queries a Differential Expression gene list for drivers in a given module biotype association
#' @description given a Differential Expression gene list it is a great interest to query the DE list for each driver identity.
#' @param DE.list a diff.expr arkas list with limma meta data
#' @param Module.color  a color module to query
#' @param path  path to the wgcnaDbLite DB
#' @param dbname  the name of the database without the .sqlite extension
#' @param module  if the dbname is null then the method expects the user to input the module to query.
#' @export
#' @return data print out to screen
queryDE_driver<-function(DE.list=NULL, Module.color="blue",path=".",dbname=NULL,module=NULL){
  stopifnot(is.null(DE.list)==FALSE) 

  if(is.null(dbname)==FALSE){
  moduleColor<-modulesBy(wgcnaDbLite(dbname),Module.color=Module.color)
  } else{
  stopifnot(is.null(module)==FALSE)
  moduleColor<-module
  }
  
  drivers<-DE.list[rownames(DE.list)%in%moduleColor$row_names,]
  up<-drivers[which(drivers$logFC>0),]
  cat(paste0("there are ",nrow(up)," +logFC DE \n"))
  down<-drivers[which(drivers$logFC<0),]
  cat(paste0("thre are ",nrow(down)," -logFC DE \n"))
  stopifnot((nrow(up)+nrow(down))==nrow(drivers))
  
  ###fix me: add a print out report
  write.csv(drivers,file=paste0(Module.color,"_DE.List.csv"),row.names=TRUE,quote=FALSE)
  return(NULL)
} #main_
