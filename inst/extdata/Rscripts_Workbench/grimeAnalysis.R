#' @title this analyzes the grime FLT3 and DNMT3a mouse model
#' @description this produces a repeat analysis of grime paper, overall we see ERVb7 family up regulated with the loss of DNMT3a which means that it is likely methylated by DNMT3a, and the loss incurs hypomethylation thus upregulating it. The MxCre is a recombinase that allows for deletion, fusion and alteration of DNMT3a. where floxed means they sandwhich the region wtih a gene with a double homozygous mutant of Flt3, and the MxCre creates a recombination mutant on one allele.  
grimeAnalysis<-function(grimes){

tt<-kexpAnalysis(grimes,byWhich="repeat",adjust="BH",read.cutoff=2)

}
