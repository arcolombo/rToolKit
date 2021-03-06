% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wgcna_scatterMod.R
\name{wgcna_scatterMod}
\alias{wgcna_scatterMod}
\title{plots a scatter correlation image of a fixed biotype and fixed eigengene module color}
\usage{
wgcna_scatterMod(lnames, biotype = c("ERV1", "ERV2", "Endogenous Retrovirus",
  "ERV3", "ERVL", "L1", "L2", "LTR Retrotransposon"), useBiCor = TRUE)
}
\arguments{
\item{lnames}{a list of the results from the call wgcna method}

\item{biotype}{the tx biotype of interest}

\item{useBiCor}{bicor is a WGCNA functoin that uses biweight midcorrelations and is robust against outliers}
}
\value{
prints images and returns module significance using student t test and fisher exact test
}
\description{
allows for exploration of data presented in wgcna_Cormap of individual modules plotted by individual traits and shows their absolute correlation.  also returns a data frame of modules' correlation and student and fisher pvalues in terms of ENSG id.  wgcna_filter is downstream of this analysis. note: you only want to call this function with plotALL=TRUE once, this will identify modules that have high absolute correlation.
}
