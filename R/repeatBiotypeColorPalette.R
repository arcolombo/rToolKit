#' @title this will uniformly handle the color palette for repeat plotting
#' @description attempts to uniformly control the color assignment for biotype and family color plots for heatmapping
#' @param kexp a kallisto experiment annotated
#' @import grDevices
#' @import colorRamps
#' @import colorspace
#' @import TxDbLite
#' @import arkas
#' @export
#' @return a vector ready for rowAnnotation ComplexHeatmap
repeatBiotypeColorPalette<-function(kexp){
kexp<-findRepeats(kexp)
color.df<-data.frame(txb=levels(factor(rowRanges(kexp)$tx_biotype)),palette="none",stringsAsFactors=FALSE)
 ##colors for acromeric/centromeric
acro.group<-length(c(grep("acromeric",color.df$txb,ignore.case=TRUE),grep("centromeric",color.df$txb  )))
color.df[  c(grep("acromeric",color.df$txb,ignore.case=TRUE),grep("centromeric",color.df$txb  ))  ,]$palette<-acro.centro(acro.group)  

 ##colors for Alu/SVA
 alu.group<-length( c(grep("Alu",color.df$txb,ignore.case=FALSE),grep("SVA",color.df$txb  )) )
 color.df[c(grep("Alu",color.df$txb,ignore.case=FALSE),grep("SVA",color.df$txb  ))  ,]$palette<-alu.sva(alu.group )

 #erv colors
 endg.group<-length(grep("Endogenous Retrovirus",color.df$txb))
 if(endg.group>0){
 color.df[c( grep("Endogenous Retrovirus",color.df$txb)),]$palette<-rgb(t(col2rgb(colors()[c(28)])),maxColorValue=255)
 }
  erv1.group<-length(which(color.df$txb=="ERV1"))
  if(erv1.group>0){
 color.df[which(color.df$txb=="ERV1"),]$palette<-rgb(t(col2rgb(colors()[c(565)])),maxColorValue=255)
  }
 erv3.group<-length(which(color.df$txb=="ERV3"))
  if(erv3.group>0){
 color.df[which(color.df$txb=="ERV3"),]$palette<-rgb(t(col2rgb(colors()[c(114)])),maxColorValue=255)
  }
 ervk.group<-length(which(color.df$txb=="ERVK"))
  if(ervk.group>0){
 color.df[which(color.df$txb=="ERVK"),]$palette<-rgb(t(col2rgb(colors()[c(553)])),maxColorValue=255)
  }
 ervl.group<-length(which(color.df$txb=="ERVL"))
  if(ervl.group>0){
   color.df[which(color.df$txb=="ERVL"),]$palette<-rgb(t(col2rgb(colors()[c(141)])),maxColorValue=255)
  }
 
 ##EUTR Colors
  eutr.group<-length(grep("EUTR",color.df$txb,ignore.case=TRUE))
  color.df[grep("EUTR",color.df$txb,ignore.case=TRUE),]$palette<-eutr(eutr.group)

  ##L1 Colors
  l1.group<-length( c(grep("L1",color.df$txb)))
  if(l1.group>0){
 color.df[ which(color.df$txb=="L1"),]$palette<-rgb(t(col2rgb(colors()[c(139)])),maxColorValue=255)
  }
 l2.group<-length( c(grep("L2",color.df$txb)))
  if(l2.group>0){
 color.df[ which(color.df$txb=="L2"),]$palette<-rgb(t(col2rgb(colors()[c(373)])),maxColorValue=255)
  }
#LTR
#Dark grey
  ltr.group<-length( c(grep("LTR Retrotransposon",color.df$txb)))
   if(ltr.group>0){
  color.df[ which(color.df$txb=="LTR Retrotransposon"),]$palette<-rgb(t(col2rgb(colors()[c(629)])),maxColorValue=255)
  }
#CR1
#magenta
  cr1.group<-length( c(grep("CR1",color.df$txb)))
   if(cr1.group>0){
  color.df[ which(color.df$txb=="CR1"),]$palette<-rgb(t(col2rgb(colors()[c(451)])),maxColorValue=255)
   }
#MULE
#bright green
   mule.group<-length( c(grep("MULE",color.df$txb)))
   if(mule.group>0){
   color.df[ which(color.df$txb=="MULE"),]$palette<-rgb(t(col2rgb(colors()[c(393)])),maxColorValue=255)
  }


  ##SAT colors
  sat.group<-length( c(grep("SAT",color.df$txb),grep("satellite",color.df$txb,ignore.case=TRUE)  )) 
  color.df[c(grep("SAT",color.df$txb),grep("satellite",color.df$txb,ignore.case=TRUE)),]$palette<-sat(sat.group)


  ##SINE
  sine.group<-length( c(grep("SINE",color.df$txb,ignore.case=FALSE),grep("snRNA",color.df$txb  ))  )
  color.df[ c(grep("SINE",color.df$txb,ignore.case=FALSE),grep("snRNA",color.df$txb  ))  ,]$palette<-sine(sine.group  )

  ##RTE RTEX
   rte.group<-length(c(grep("Repetitive element",color.df$txb),grep("RTE",color.df$txb,ignore.case=FALSE)  )  )
  color.df[c(grep("Repetitive element",color.df$txb),grep("RTE",color.df$txb,ignore.case=FALSE)  ),]$palette<-rte.re(rte.group)


 ##Transposable elements
 #brown/green and light grey
    trans.group<-length( c(grep("DNA transposon",color.df$txb),grep("Transposable Element",color.df$txb)  ))
   color.df[c(grep("DNA transposon",color.df$txb),grep("Transposable Element",color.df$txb)  )  ,]$palette<-trans(trans.group  )

 ##bright green and dark purple
    tc.group<-length( c(grep("TcMar",color.df$txb),grep("Mariner/Tc1",color.df$txb)  ))
   if(tc.group>0){
  color.df[c(grep("Mariner/Tc1",color.df$txb),grep("TcMar",color.df$txb)),]$palette<-rgb(t(col2rgb(colors()[c(254,551)])),maxColorValue=255)
  }
 
 ##red and orange
 mer.group<-length( c(grep("Merlin",color.df$txb),grep("MIR",color.df$txb)  ))
  if(mer.group==1){
   color.df[c(grep("Merlin",color.df$txb),grep("MIR",color.df$txb)),]$palette<-rgb(t(col2rgb(colors()[c(503)])),maxColorValue=255)
  }else if(mer.group==2){
   color.df[c(grep("Merlin",color.df$txb),grep("MIR",color.df$txb)),]$palette<-rgb(t(col2rgb(colors()[c(503,34)])),maxColorValue=255)
 }
 

  ##remnants
 remnant.group<-length(which(color.df$palette=="none")  )
  if(remnant.group>0){
  color.df[which(color.df$palette=="none"),]$palette<-rainbow(remnant.group )
  }
 tx.color.master<-color.df$palette
names(tx.color.master)<-color.df$txb
# txb.sel.col=list(transcript_biotype=tx.color.master)
#id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype
#rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))
 return(tx.color.master)

}

#' @import grDevices
#' @import colorRamps
#' @import colorspace
acro.centro<-function (n, h = 360, c. = c(53, 100), l = c(97, 86), power = 3, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}




#' @import grDevices
#' @import colorRamps
#' @import colorspace
alu.sva<-function (n, h = 245, c. = c(63, 89), l = c(83, 100), power = 3, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}


#' @import grDevices
#' @import colorRamps
#' @import colorspace
ervs<-function (n, h = 260, c. = c(100, 10), l = c(30, 100), power = 0.7, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}


#' @import grDevices
#' @import colorRamps
#' @import colorspace
eutr<-function (n, h = -275, c. = c(100, 0), l = c(22, 83), power = 0.840909090909091, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}

#' @import grDevices
#' @import colorRamps
#' @import colorspace
l1<-function (n, h = 295, c. = c(89, 3), l = c(38, 54), power = 0.818181818181818, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}

#' @import grDevices
#' @import colorRamps
#' @import colorspace
sat<-function (n, h = -194, c. = c(100, 71), l = c(21, 25), power = 0.159090909090909, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}

#' @import grDevices
#' @import colorRamps
#' @import colorspace
sine<-function (n, h = 33, c. = c(80, 10), l = c(30, 83), power = 0.7, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}


#' @import grDevices
#' @import colorRamps
#' @import colorspace
rte.re<-function (n, h = 3, c. = c(0, 0), l = c(2, 47), power = 1.32954545454545, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}

#' @import grDevices
#' @import colorRamps
#' @import colorspace
trans<-function (n, h = 68, c. = c(39, 10), l = c(51, 79), power = 0.4, 
    fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
        ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}



