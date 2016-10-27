#' @title use the probability of a fit distribution in simulated data for exp dist
#' @description use the probability of a single observation as the true population parameter and simulates logTPM values based on a fitted exp distribution as the population dist.  takes the observation probabilty and adjusts by an epsilon to then remap onto a tpm value.
#' @param dum  the simulated data frame
#' @param delta.max the integer to adjust the probability values by
#' @param lamb  lambda rate
#' @export
#' @return a data frame of logTPM values
simCorrect<-function(dum=NULL,delta.max=0.05,lamb=fit.exp$estimate["rate"]){
stopifnot(is.null(dum)==FALSE)
DC0h_mm1<-dum$DC0h_mm1
DC0h_mm05.prob<-dum$true.prob+delta.max/2
if(any(DC0h_mm05.prob>1)){
DC0h_mm05.prob[which(DC0h_mm05.prob>1)]<-dum$true.prob[which(DC0h_mm05.prob>1)]
}
DC0h_mm05.logTpm<-qexp(DC0h_mm05.prob,rate=lamb)

DC0h_mmN05.prob<-dum$true.prob-delta.max
if(any(DC0h_mmN05.prob<0)){
DC0h_mmN05.prob[which(DC0h_mmN05.prob<0)]<-dum$true.prob[which(DC0h_mmN05.prob<0)]
}
DC0h_mmN05.logTpm<-qexp(DC0h_mmN05.prob,rate=lamb)


dum<-data.frame(DC0h_mm1=DC0h_mm1,
                  DC0h_mm05=DC0h_mm05.logTpm,
                  DC0h_mmN05=DC0h_mmN05.logTpm,
                  true.prob=dum$true.prob,
                  DC0h_mm05.prob=DC0h_mm05.prob,
                  DC0h_mmN05.prob=DC0h_mmN05.prob,
                  delta=DC0h_mm05.prob-dum$true.prob,
                  delta2=DC0h_mmN05.prob-dum$true.prob)


####



return(dum)

}#main
