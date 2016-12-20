#' @title this will uniformly handle the color palette for repeat plotting
#' @description attempts to uniformly control the color assignment for biotype and family color plots for heatmapping
#' @param kexp a kallisto experiment annotated
#' @import grDevices
#' @import colorRamps
#' @export
#' @return a vector ready for rowAnnotation ComplexHeatmap
repeatColorPalette<-function(kexp){

color.df<-data.frame(txb=levels(factor(rowRanges(kexp)$tx_biotype)),palette="none",stringsAsFactors=FALSE)

color.df[color.df$txb%in%c("acromeric","centromeric"),]$palette<-acro.centro(2)
 color.df[color.df$txb%in%c("Alu","SVA"),]$palette<-alu.sva(2)
 color.df[6:10,]$palette<-ervs(5)
 color.df[11:37,]$palette<-eutr(27)
color.df[c(43,44,45),]$palette<-l1(3)
color.df[57:59,]$palette<-sat(3)
color.df[c(60:61,62),]$palette<-sine(3)
color.df[c(54,55,56,51),]$palette<-rte.re(4)
 color.df[c(5,66),]$palette<-trans(2)

color.df[which(color.df$palette=="none"),]$palette<-rainbow(17)

 tx.color.master<-color.df$palette
names(tx.color.master)<-color.df$txb
# txb.sel.col=list(transcript_biotype=tx.color.master)
#id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype
#rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))
 return(tx.color.master)

}

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



