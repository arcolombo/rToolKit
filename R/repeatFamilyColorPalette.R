#' @title this will uniformly handle the color palette for repeat plotting
#' @description attempts to uniformly control the color assignment for biotype and family color plots for heatmapping
#' @param kexp a kallisto experiment annotated
#' @import grDevices
#' @import colorRamps
#' @importFrom hex colorspace
#' @import TxDbLite
#' @import arkas
#' @export
#' @return a vector ready for rowAnnotation ComplexHeatmap
repeatFamilyColorPalette<-function(kexp){

  kexp<-findRepeats(kexp)
  color.df<-data.frame(txb=levels(factor(rowRanges(kexp)$gene_biotype)),palette="none",stringsAsFactors=FALSE)


 ##colors for dna element
  dna.group<-length(grep("DNA_element",color.df$txb,ignore.case=TRUE))
  color.df[  c(grep("DNA_element",color.df$txb,ignore.case=TRUE)),]$palette<-dna.element(dna.group ) ##yellow  

 ##colors for line
 line.group<-length( c(grep("LINE",color.df$txb,ignore.case=FALSE)) )
 color.df[c(grep("LINE",color.df$txb,ignore.case=FALSE))  ,]$palette<-line.family(line.group )

 #erv colors
 erv.group<-length(c(grep("ERV",color.df$txb),grep("Endogenous Retrovirus",color.df$txb)) )
 color.df[c( grep("Endogenous Retrovirus",color.df$txb),grep("ERV",color.df$txb)),]$palette<-ervs(erv.group)

 ##LTR Colors
  ltr.group<-length(grep("LTR_element",color.df$txb))
  color.df[grep("LTR_element",color.df$txb),]$palette<-ltr.family(ltr.group)

  ##SINE
  sine.group<-length( c(grep("SINE",color.df$txb,ignore.case=FALSE))  )
color.df[ c(grep("SINE",color.df$txb,ignore.case=FALSE) )  ,]$palette<-sine.family(sine.group  )



  ##remnants
 remnant.group<-length(which(color.df$palette=="none")  )
color.df[which(color.df$palette=="none"),]$palette<-other.family(remnant.group )

 gene.color.master<-color.df$palette
names(gene.color.master)<-color.df$txb
# txb.sel.col=list(transcript_biotype=tx.color.master)
#id2<-names(txb.sel.col[[1]])%in%rowRanges(kexp)[rownames(rpt.mt)]$tx_biotype
#rA<-rowAnnotation(txb.df,col=list(transcript_biotype=(txb.sel.col[[1]][id2])))
 return(gene.color.master)

}


dna.element<-function (n, h = -327, c. = c(100, 100), l = c(72, 87), power = 0.568181818181818, 
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


line.family<-function (n, h = 280, c. = c(100, 95), l = c(50, 100), power = 0, 
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
 

ltr.family<-function (n, h = 0, c. = c(0, 0), l = c(0, 0), power = 0, fixup = TRUE, 
    gamma = NULL, alpha = 1, ...) 
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



sine.family<-function (n, h = 245, c. = c(100, 12), l = c(41, 57), power = 0, 
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



other.family<-function (n, h = 135, c. = c(100, 100), l = c(64, 100), power = 1.27272727272727, 
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

