#' @title this analyzes hubs that are preserved
kmeWiseAnalysis<-function(gnames,rnames){


#FIX ME: the problem here is that repeats and genes do not share any common genes so a consensus I cant figure out how to do.
 ###One method is to take the signedKME(gEXP,gMEs) and plot against signedKME(gEXP,repeatMEs)   this will take the gene connectivity plotted against how well the genes are connected to repeat eigenValues.   not sure yet how to do this.
###hub of genes and repeats.
  kME<-signedKME(datExpr,MEs)
  MMPvalue<-corPvalueStudent(as.matrix(kME),nrow(datExpr)) ##nSampes=nrow
  name.id<-match(rownames(kME),rownames(MMPvalue))
  mmpValue<-as.data.frame(MMPvalue[name.id,])
  stopifnot(all(rownames(kME)==rownames(mmpValue))==TRUE)
###hub of genes and repeats.
  kME<-signedKME(datExpr,MEs)
  MMPvalue<-corPvalueStudent(as.matrix(kME),nrow(datExpr)) ##nSampes=nrow
  name.id<-match(rownames(kME),rownames(MMPvalue))
  mmpValue<-as.data.frame(MMPvalue[name.id,])
  stopifnot(all(rownames(kME)==rownames(mmpValue))==TRUE)
   kmeRed<-data.frame(row_names=rownames(kME),
                      kMEred=kME$kMEred,
                      p.Value=mmpValue$kMEred)
   redSig<-kmeRed[which(kmeRed$p.Value<0.05),]
   redSig2<-redSig[which(redSig$kMEred>0.7),]

 ###now for repeats

 ###hub of genes and repeats.
  rkME<-signedKME(rdatExpr,rMEs)
  rMMPvalue<-corPvalueStudent(as.matrix(rkME),nrow(rdatExpr)) ##nSampes=nrow
  rname.id<-match(rownames(rkME),rownames(rMMPvalue))
  rmmpValue<-as.data.frame(rMMPvalue[rname.id,])
  stopifnot(all(rownames(rkME)==rownames(rmmpValue))==TRUE)
   rkmeBlue<-data.frame(row_names=rownames(rkME),
                      kMEblue=rkME$kMEblue,
                      p.Value=rmmpValue$kMEblue)
   blueSig<-rkmeBlue[which(rkmeBlue$p.Value<0.05),]
   blueSig2<-blueSig[which(blueSig$kMEred>0.7),]


   kmeRed<-data.frame(row_names=rownames(kME),
                      kMEred=kME$kMEred,
                      p.Value=mmpValue$kMEred)
   redSig<-kmeRed[which(kmeRed$p.Value<0.05),]
   redSig2<-redSig[which(redSig$kMEred>0.7),]

 ###now for repeats

 ###hub of genes and repeats.
  rkME<-signedKME(rdatExpr,rMEs)
  rMMPvalue<-corPvalueStudent(as.matrix(rkME),nrow(rdatExpr)) ##nSampes=nrow
  rname.id<-match(rownames(rkME),rownames(rMMPvalue))
  rmmpValue<-as.data.frame(rMMPvalue[rname.id,])
  stopifnot(all(rownames(rkME)==rownames(rmmpValue))==TRUE)
   rkmeBlue<-data.frame(row_names=rownames(rkME),
                      kMEblue=rkME$kMEblue,
                      p.Value=rmmpValue$kMEblue)
   blueSig<-rkmeBlue[which(rkmeBlue$p.Value<0.05),]
   blueSig2<-blueSig[which(blueSig$kMEred>0.7),]


  ###correlate top ranked 



 ##plot softConnectivity_genes vs softConnectivity_repeats



}
