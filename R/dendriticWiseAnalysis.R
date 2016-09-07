#' @title analyzing the dendritic tet2 KO experiment requires simulating replicates
#' @description the tet2 KO has no replicates so this method will simulate samples using a uniform distribution to create multiple replicates with some technical digital variation.  a DE call will be used but the results are from simulated data.
#' @param kexp a kallisto experiment repeat stage
#' @param bundleID tx_id or gene_id
#' @param read.cutoff floor threshold
#' @param whichModel norm pois gamma exp or unif, these are options to test the fit model against
#' @param n.sims  5e4 numeric to simulate trials total
#' @param minTPM for exp model TPMs are used and this is the floor TPM cut
#' @param newCols integer, new columns to simulate: breathe in the newness of life.
#' @import fitdistrplust
#' @import logspline
#' @export
#' @return a matrix/df of simulated columns dim n X (2*newCols+2)
dendriticWiseAnalysis<-function(kexp,bundleID="tx_id",read.cutoff=1,whichModel=c("norm","pois","gamma","exp","unif"),n.sims=5e4,minTPM=0.3,newCols=3){
  ##FIX ME : add uniform simulations of data
   Oh<-grep("DC0h_mm",colnames(kexp))
   Fourh<-grep("DC4h_mm",colnames(kexp)) 
   kexp<-findRepeats(kexp)
   bundledCounts<-collapseBundles(kexp,bundleID=bundleID,read.cutoff=read.cutoff)
   x<-(bundledCounts[,Oh])
   y<-(bundledCounts[,Fourh])

  fit.norm<-fitdist((x),"norm")
  fit.poi<-fitdist(round(x),"pois")
  fit.gamma<-fitdist(x,"gamma",method="mme")
  message("norm")
  plot(fit.norm)
  readkey()
  message("poi")
  plot(fit.poi)
  readkey()
  message("gamma")
  plot(fit.gamma)
  readkey()


  if(whichModel=="exp"){
  message("boostrapping and simulations are done on log 1+tpm_ij")
  bundledCounts<-collapseTpm(kexp,bundleID=bundleID,minTPM=minTPM)
  x<-bundledCounts[,Oh]
  y<-bundledCounts[,Fourh]
  x<-log(1+x)
  x<-x[which(x>0)]
  fit.exp<-fitdist(x,"exp")
  message('exp')
  plot(fit.exp)
  readkey()
  
stats <- replicate(n.sims, {
  r <- rexp(n = length(x),
                 rate= fit.exp$estimate["rate"] )
  as.numeric(ks.test(r,
                    "pexp",
                     rate= fit.exp$estimate["rate"] ))
})
  plot(ecdf(stats), las = 1, main = "KS-test statistic simulation (CDF)", col = "darkorange", lwd = 1.7)
  grid()
  print(fit.exp$aic)
  fit <- logspline(stats)  ###many NAs
##bootstrapping
  xs <- seq(min(x),max(x),len=length(x))
  true.exp <- rexp(1e6, rate= fit.exp$estimate["rate"])
  boot.pdf <- sapply(1:1000, function(i) {
  xi <- sample(x, size=length(x), replace=TRUE)
  MLE.est <- fitdist(xi, distr="exp")
  dexp(xs, rate=MLE.est$estimate["rate"])})

boot.cdf <- sapply(1:1000, function(i) {
  xi <- sample(x, size=length(x), replace=TRUE)
  MLE.est <- suppressWarnings(fitdist(xi, distr="exp"))
  pexp(xs, rate= MLE.est$estimate["rate"]  )})
#-----------------------------------------------------------------------------
# Plot PDF
#-----------------------------------------------------------------------------
  par(bg="white", las=1, cex=1.2)
  plot(xs, boot.pdf[, 1], type="l", col=rgb(.6, .6, .6, .1), ylim=range(boot.pdf),
     xlab="x", ylab="Probability density")
for(i in 2:ncol(boot.pdf)) lines(xs, boot.pdf[, i], col=rgb(.6, .6, .6, .1))

# Add pointwise confidence bands
  quants <- apply(boot.pdf, 1, quantile, c(0.025, 0.5, 0.975))
  min.point <- apply(boot.pdf, 1, min, na.rm=TRUE)
  max.point <- apply(boot.pdf, 1, max, na.rm=TRUE)
  lines(xs, quants[1, ], col="red", lwd=1.5, lty=2)
  lines(xs, quants[3, ], col="red", lwd=1.5, lty=2)
  lines(xs,quants[2,],col="darkred",lwd=2)
#----------------------------------------------------------------------------
  par(bg="white", las=1, cex=1.2)
  plot(xs, boot.cdf[, 1], type="l", col=rgb(.6, .6, .6, .1), ylim=range(boot.cdf),
     xlab="x", ylab="F(x)")
for(i in 2:ncol(boot.cdf)) lines(xs, boot.cdf[, i], col=rgb(.6, .6, .6, .1))

# Add pointwise confidence bands
  quants <- apply(boot.cdf, 1, quantile, c(0.025, 0.5, 0.975))
  min.point <- apply(boot.cdf, 1, min, na.rm=TRUE)
  max.point <- apply(boot.cdf, 1, max, na.rm=TRUE)
  lines(xs, quants[1, ], col="red", lwd=1.5, lty=2)
  lines(xs, quants[3, ], col="red", lwd=1.5, lty=2)
  lines(xs, quants[2, ], col="darkred", lwd=2)
} else if(whichModel=="norm"){

  stats <- replicate(n.sims, {      
  r <- rnorm(n = length(x),
                 mean= fit.norm$sd["mean"],
                 sd = fit.norm$sd["sd"] )
  as.numeric(ks.test(r,
                    "pnorm",
                     mean= fit.norm$sd["mean"],
                     sd = fit.norm$sd["sd"])  )      
})
  plot(ecdf(stats), las = 1, main = "KS-test statistic simulation (CDF)", col = "darkorange", lwd = 1.7)
  grid()
   ##pvalue using simulated distribution
  fit <- logspline(stats)
1 - plogspline(ks.test(x,
                       "pnorm",
                       mean= fit.norm$sd["mean"],
                      sd = fit.norm$sd["sd"]),
                      fit)
#################
 stats <- replicate(n.sims, {
  r <- rpois(n = length(x),
                 lambda= fit.poi$estimate["lambda"] )
  as.numeric(ks.test(r,
                    "ppois",
                     lambda= fit.poi$estimate["lambda"] ))
})
 plot(ecdf(stats), las = 1, main = "KS-test statistic simulation (CDF)", col = "darkorange", lwd = 1.7)
grid()

  ##pvalue using simulated distribution
  fit <- logspline(stats)  ###many NAs
##bootstrapping
  xs <- seq.int(65, 500,1)

  true.pois <- rpois(1e6, lambda= fit.poi$estimate["lambda"])
  boot.pdf <- sapply(1:1000, function(i) {
  xi <- sample(round(x), size=length(x), replace=TRUE)
  MLE.est <- fitdist(xi, distr="pois")
  dpois(xs, lambda=MLE.est$estimate["lambda"])})

boot.cdf <- sapply(1:1000, function(i) {
  xi <- sample(round(x), size=length(x), replace=TRUE)
  MLE.est <- suppressWarnings(fitdist(xi, distr="pois"))  
  ppois(xs, lambda= MLE.est$estimate["lambda"]  )})   

#-----------------------------------------------------------------------------
# Plot PDF
#-----------------------------------------------------------------------------

par(bg="white", las=1, cex=1.2)
plot(xs, boot.pdf[, 1], type="l", col=rgb(.6, .6, .6, .1), ylim=range(boot.pdf),
     xlab="x", ylab="Probability density")
for(i in 2:ncol(boot.pdf)) lines(xs, boot.pdf[, i], col=rgb(.6, .6, .6, .1))

# Add pointwise confidence bands

quants <- apply(boot.pdf, 1, quantile, c(0.025, 0.5, 0.975))
min.point <- apply(boot.pdf, 1, min, na.rm=TRUE)
max.point <- apply(boot.pdf, 1, max, na.rm=TRUE)
lines(xs, quants[1, ], col="red", lwd=1.5, lty=2)
lines(xs, quants[3, ], col="red", lwd=1.5, lty=2)

DC0h_mm10<-vector()
DC0h_mm15<-vector()
DC0h_mm20<-vector()

DC0h_mm1<-x  ##log TPM of original data
#DC0h_mm2<-rexp(length(x),rate=fit.exp$estimate["rate"]) ##simulated under sample distr
true.prob<-pexp(x,rate=fit.exp$estimate["rate"])
#sample.prob<-pexp(DC0h_mm2,rate=fit.exp$estimate["rate"])
simDF<-data.frame(DC0h_mm1=DC0h_mm1,
              #    DC0h_mm2=DC0h_mm2,
                  true.prob=true.prob)
               #   sample.prob=sample.prob,
               #   delta=sample.prob-true.prob )

simCol_05<-simCorrect(simDF,delta.max=0.103,delta.min=-0.103,lamb=as.numeric(fit.exp$estimate["rate"]),x=DC0h_mm1,true.prob=true.prob)
simCol_07<-simCorrect(simDF,delta.max=0.106,delta.min=-0.106,lamb=as.numeric(fit.exp$estimate["rate"]),x=DC0h_mm1,true.prob=true.prob)
simCol_08<-simCorrect(simDF,delta.max=0.109,delta.min=-0.109,lamb=as.numeric(fit.exp$estimate["rate"]),x=DC0h_mm1,true.prob=true.prob)
simCol_10<-simCorrect(simDF,delta.max=0.11,delta.min=-0.11,lamb=as.numeric(fit.exp$estimate["rate"]),x=DC0h_mm1,true.prob=true.prob)
simCol_15<-simCorrect(simDF,delta.max=0.15,delta.min=-0.15,lamb=as.numeric(fit.exp$estimate["rate"]),x=DC0h_mm1,true.prob=true.prob)

stopifnot(all(x==simCol_05$DC0h_mm1))
stopifnot(all(x==simCol_10$DC0h_mm1))
stopifnot(all(x==simCol_15$DC0h_mm1))


DC0h_mm<-data.frame(DC0h_mm1=x,
                    DC0h_mm05=simCol_05$DC0h_mm2,
                    DC0h_mm10=simCol_10$DC0h_mm2,
                    DC0h_mm07=simCol_07$DC0h_mm2,
                    DC0h_mm08=simCol_08$DC0h_mm2,
                    DC0h_mm15=simCol_15$DC0h_mm2,stringsAsFactors=FALSE)


#-----------------------------------------------------------------------------
  par(bg="white", las=1, cex=1.2)
  plot(xs, boot.cdf[, 1], type="l", col=rgb(.6, .6, .6, .1), ylim=range(boot.cdf),
     xlab="x", ylab="F(x)")
for(i in 2:ncol(boot.cdf)) lines(xs, boot.cdf[, i], col=rgb(.6, .6, .6, .1))

# Add pointwise confidence bands

quants <- apply(boot.cdf, 1, quantile, c(0.025, 0.5, 0.975))
min.point <- apply(boot.cdf, 1, min, na.rm=TRUE)
max.point <- apply(boot.cdf, 1, max, na.rm=TRUE)
lines(xs, quants[1, ], col="red", lwd=1.5, lty=2)
lines(xs, quants[3, ], col="red", lwd=1.5, lty=2)
lines(xs, quants[2, ], col="darkred", lwd=2)

} else if(whichModel=="unif"){
 ##sample each row uniformly to regenerate the data uniformly.
   Oh<-grep("DC0h_mm",colnames(kexp))
   Fourh<-grep("DC4h_mm",colnames(kexp))
   kexp<-findRepeats(kexp)
   bundledCounts<-collapseBundles(kexp,bundleID=bundleID,read.cutoff=read.cutoff)
   x<-(bundledCounts[,Oh])
   y<-(bundledCounts[,Fourh])
  qbin<-quantile(x)  


} else {
 message("currently support exp and norm")
}



} ##{main
