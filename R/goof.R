#' @title goodness of fit for repeat data to poisson model or linear
#' @description goodness of fit for repeat data follow simulating from a poisson and from normal distribution.  It first looks at all the samples together as one observation, then examines each sample individually to see how each sample fits.  the concern is that if all the samples are examined 
#' @param kexp a kallisto experiment of repeat stage or all transcript stage
#' @param n the number of samples to draw from the poisson distribution
#' @param numberOut the number to plot out
#' @param n.trials the number of trials to create and sample from
#' @param datamat a matrix like object in the case where a kexp is undesired
#' @export
#' @return images to test poisson hyp
goof<-function(kexp,n=100,numberOut=999,n.trials=10000,dataMat=NULL){

###FIX ME: need to test GOF to the correct model.  how well does the repeat data fit to a poisson?  how well does the repeat data fit to a Guassian model?  look at poissonSeq by li et. al

##notes: as n, the sample size increases the data will approximate a normal distribution. 
##Simulates poisson dist from the observed sample statistics
  if(is.null(dataMat)==TRUE){
  mu<-mean(counts(kexp))
  } else if(is.null(dataMat)==FALSE){
  mu<-mean(dataMat)
  } 
  n<-n
  n.trials<-n.trials
  sim<-replicate(n.trials,{x<-rpois(n,mu);c(mean(x),var(x))})
  xy<-apply(sim,1,function(x) quantile(x,probs=seq(0.001,0.999,length.out=numberOut)))
  plot(xy,xlab="mean",ylab="var",col=hsv(0,0,0.25,.4),main=paste0("Poisson Random Generation from Repeat mean ",mu," s.size:",n," trials:",n.trials))
  readkey()
  for(i in 1:ncol(counts(kexp))){
  if(is.null(dataMat)==TRUE){ 
  sampleX<-counts(kexp)[,i]
  } else if(is.null(dataMat)==FALSE){
   sampleX<-dataMat[,i]
   }
   mu<-mean(sampleX)
    n<-n
  n.trials<-n.trials
  sim<-replicate(n.trials,{x<-rpois(n,mu);c(mean(x),var(x))})
  xy<-apply(sim,1,function(x) quantile(x,probs=seq(0.001,0.999,length.out=numberOut)))
  plot(xy,xlab="mean",ylab="var",col=hsv(0,0,0.25,.4),main=paste0("Simulated Poisson from ",colnames(kexp)[i]," mean: ",mu," s.size:",n," trials:",n.trials ))
  readkey()
  }

#####FIX ME ::: add chi^2 fit  
###fix me:: add a linear model regression 
 
###### poisson seems to have a goodfit call a method to find goodfit stat.... examine stat......  if pass, then PoissonSeq, and cqnDE are allowed 




}
