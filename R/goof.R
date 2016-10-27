#' @title goodness of fit for repeat data to poisson model or linear
#' @description goodness of fit for repeat data follow simulating from a poisson and from normal distribution.  It first looks at all the samples together as one observation, then examines each sample individually to see how each sample fits.  the concern is that if all the samples are examined 
#' @param kexp a kallisto experiment of repeat stage or all transcript stage
#' @param n the number of samples to draw from the poisson distribution
#' @param numberOut the number to plot out
#' @param n.trials the number of trials to create and sample from
#' @param datamat a matrix like object in the case where a kexp is undesired
#' @param bySamples boolean if true the good fit is calculate 
#' @import limma
#' @import vcd
#' @import tseries
#' @import nortest
#' @export
#' @return images to test poisson hyp
goof<-function(kexp,n=100,numberOut=999,n.trials=10000,dataMat=NULL,bySamples=FALSE){
##notes: as n, the sample size increases the data will approximate a normal distribution. 

##Simulates poisson dist from the observed sample statistics
  if(is.null(dataMat)==TRUE){
  mu<-mean(counts(kexp))
  cnts<-counts(kexp)
  } else if(is.null(dataMat)==FALSE){
  mu<-mean(dataMat)
  cnts<-dataMat
  } 
  n<-n
  n.trials<-n.trials
  sim<-replicate(n.trials,{x<-rpois(n,mu);c(mean(x),var(x))})
  xy<-apply(sim,1,function(x) quantile(x,probs=seq(0.001,0.999,length.out=numberOut)))
  plot(xy,xlab="mean",ylab="var",col=hsv(0,0,0.25,.4),main=paste0("Poisson Random Generation from Repeat mean ",mu," s.size:",n," trials:",n.trials))
 
  z<-xy[,1]^2
  fit<-lm(xy[,2]~xy[,1]+z)
  x0<-seq(mu-2*sqrt(mu/n),mu+2*sqrt(mu/n),length.out=99)
  u<-cbind(rep(1,length(x0)),x0,x0^2)%*%coefficients(fit)
  lines(x0,u,lwd=2,col=hsv(0,.8,.8,1))
  message(paste0("quadratic coefficient ",coefficients(fit)[3])) 
  q <- (1:100 - 0.5) / 100
  lines(qnorm(q,mean=mu,sd=sqrt(mu/n)),mu*qchisq(q,df=n-1)/(n-1),col="Blue")
  readkey() 

  y<-rnbinom(n,mu,prob=0.5)
  x<-rpois(n,mu)
  fit.y<-glm(y~1,family=poisson())
  fit.x<-glm(x~1,family=poisson())
  message("plotting neg.binom")
  plot(fit.y)
  message("plotting poisson")
  plot(fit.x) 
   readkey() 
 
  x.poi<-rpois(n=n,lambda=mean(cnts))
  lambda.est<-mean(x.poi)
  tab.os<-table(x.poi)

  freq.os<-vector()
  for(i in 1: length(tab.os)) {
   freq.os[i]<-tab.os[[i]]  ## vector of emprical frequencies
   }
#  freq.ex<-(dpois(0:max(x.poi),lambda=lambda.est))*200 ## vector of fitted (expected) frequencies

 # acc <- mean(abs(freq.os-trunc(freq.ex))) ## absolute goodness of fit index acc
 # gofIndex<-acc/mean(freq.os)*100 ## relative (percent) goodness of fit index
 # message(paste0("the goodness of fit index from estimating the population parameter using sample statistic for Poisson: ",gofIndex))
  h <- hist(x.poi ,breaks=15)
  xhist <- c(min(h$breaks),h$breaks)
  yhist <- c(0,h$density,0)
  xfit <- seq(min(x.poi),max(x.poi),length=40)
  yfit <- dpois(xfit,lambda=lambda.est)
  plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), main="Poison density and histogram")
  lines(xfit,yfit, col="red")

#Perform the chi-square goodness of fit test 
  gf <- goodfit(x.poi,type= "poisson",method= "MinChisq")
  summary(gf)
  plot(gf,main="Count data vs Poisson distribution")

######testing for normality on sample distribution
  print(jarque.bera.test(x.poi))
  readkey()
  print(shapiro.test(x.poi))
  readkey()
  if(n>=7){
  print(ad.test(x.poi))
  readkey()
  print(cvm.test(x.poi))
  readkey() 
   }
  if(n>4){
  print(lillie.test(x.poi))
  readkey()
  }
  print(pearson.test(x.poi))   
  readkey()

  ###testing for normality on raw data
  if(bySamples==TRUE){
  for(i in 1:ncol(cnts)){
   x.poi<-rpois(n=n,lambda=mean(cnts[,i]))
   print(jarque.bera.test(cnts[,i]))
   readkey()
   print(shapiro.test(cnts[,i]))
   readkey()
   print(ad.test(cnts[,i]))
   readkey()
   print(cvm.test(cnts[,i]))
   readkey()
   print(lillie.test(cnts[,i]))
   readkey()
   print(pearson.test(cnts[,i]))
   readkey()

  }
 }

 if(bySamples==TRUE){
  for(i in 1:ncol(cnts)){
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

}

#####FIX ME ::: add chi^2 fit  
###fix me:: add a linear model regression 
 
###### poisson seems to have a goodfit call a method to find goodfit stat.... examine stat......  if pass, then PoissonSeq, and cqnDE are allowed 




} ##main
