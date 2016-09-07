#' @title corrects the probability differences in simulated data for exp dist
#' @description read title
simCorrect<-function(simDF=dum,delta.max=0.05,delta.min=-0.05,lamb=fit.exp$estimate["rate"],x=DC0h_mm1,true.prob=true.prob){

DC0h_mm2<-rexp(length(x),rate=lamb) ##simulated under sample distr
sample.prob<-pexp(DC0h_mm2,rate=lamb)
print(head(DC0h_mm2))
dum<-data.frame(DC0h_mm1=DC0h_mm1,
                  DC0h_mm2=DC0h_mm2,
                  true.prob=true.prob,
                  sample.prob=sample.prob,
                  delta=sample.prob-true.prob )


i<-1
while(any(dum$delta>delta.max)==TRUE){
 print(i)
 id<-which(dum$delta>delta.max)
 dum$sample.prob[id]<-pexp((dum$DC0h_mm2[id]- qexp(dum$delta[id],rate=lamb) ),rate=lamb)
 dum$DC0h_mm2[id]<-qexp(dum$sample.prob[id],rate=lamb)
 dum$delta<-dum$sample.prob-dum$true.prob
 i<-i+1
} #while 1
####
i<-1
while(any(dum$delta<(delta.min))==TRUE){
 print(i)
id<-which(dum$delta<(delta.min))
 num<-1+dum$delta[id]
 dum$sample.prob[id]<-pexp(dum$DC0h_mm2[id]+ qexp(num,rate=lamb),rate=lamb)
 dum$DC0h_mm2[id]<-qexp(dum$sample.prob[id],rate=lamb)
 dum$delta<-dum$sample.prob-dum$true.prob
 i<-i+1
} ##while 2

message(paste0("any greater ",delta.max," ",any(dum$delta)>delta.max))
message(paste0("any less ",delta.min," ",any(dum$delta)<delta.min))



return(dum)

}#main
