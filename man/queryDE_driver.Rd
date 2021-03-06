% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queryDE_driver.R
\name{queryDE_driver}
\alias{queryDE_driver}
\title{queries a Differential Expression gene list for drivers in a given module biotype association}
\usage{
queryDE_driver(DE.list = NULL, Module.color = "blue", path = ".",
  dbname = NULL, module = NULL)
}
\arguments{
\item{DE.list}{a diff.expr arkas list with limma meta data}

\item{Module.color}{a color module to query}

\item{path}{path to the wgcnaDbLite DB}

\item{dbname}{the name of the database without the .sqlite extension}

\item{module}{if the dbname is null then the method expects the user to input the module to query.}
}
\value{
data print out to screen
}
\description{
given a Differential Expression gene list it is a great interest to query the DE list for each driver identity.
}
