require(MASS)

###need to simulate the raw correlation values that should have a controlled covariance with the ranks?  
##need to simulate the ranks 
 
 a<-rnorm(n=2)
 b<-rnorm(n=2)
# c<-rnorm(n=5)
# d<-rnorm(n=5)
# e<-rnorm(n=5)
##let X be the random variable for correlations
# X<-cbind(a,b,c,d,e)
 X<-cbind(a,b)
 a1<-rnorm(n=2,mean=2,sd=1)
 b1<-rnorm(n=2,mean=2,sd=1)
# c1<-rnorm(n=5)
# d1<-rnorm(n=5)
# e1<-rnorm(n=5)
##let Y be the random variable for ranks
# Y<-cbind(a1,b1,c1,d1,e1)
  Y<-cbind(a1,b1)
  weight<-pnorm(log(Y),mean=mean(log(Y)),sd=sd(log(Y)),lower.tail=TRUE)

 ## let Z be the random bivariable for ranked correlations
Z = X*weight
#### Z is the true value

muz<-matrix(nrow=ncol(Z))
for(i in 1:ncol(Z)){
muz[i]<-mean(Z[,i])
}
stdz<-matrix(nrow=ncol(Z))
for(i in 1:ncol(Z)){
stdz[i]<-sd(Z[,i])
}
rownames(stdz)<-colnames(Z)

stats<-cbind(muz,stdz)
colnames(stats)<-c("mean","std")
rownames(muz)<-colnames(Z)

#######now to sample from theoretical distributions


###sample from the distribution, and not from the matrix Z 
require(fitdistrplus)
s<-fitdist(as.vector(Z[,1]),"norm")
s2<-fitdist(as.vector(Z[,2]),"norm")
z_sim<-rnorm(n=length(Z[,1]),mean=s$estimate["mean"],sd=s$estimate["sd"])
z_sim2<-rnorm(n=length(Z[,2]),mean=s2$estimate["mean"],sd=s2$estimate["sd"])

simulated_Z<-data.frame(z_sim,z_sim2)
colnames(simulated_Z)<-colnames(Z)
##for each vector fit a distribution



#### simulate X and simulate weights Y 
##the algorithm adjusts Y to ensure normality of X*Y

###simulate X


###simulate Y 


##adjust Y



###convolve some noise



##the task is to validate Z from simulated X and simulated Y convolved with noise.

######

