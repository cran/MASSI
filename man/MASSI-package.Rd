\name{MASSI-package}
\alias{MASSI-package}
\alias{MASSI}
\docType{package}
\title{MASSI: MicroArray Sample Sex Identifier}
\description{
MASSI uses the Y chromosome probe information to cluster and assign sex to each sample in gene expression microarray datasets.
}
\details{
\tabular{ll}{
Package: \tab MASSI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-06-11\cr
License: \tab GPL-3\cr
}

MASSI requires a sample/probe expression matrix and a list of probe id's that correspond to Y chromosome genes. 
The "massi.check" function will extract the Y chromosome probe information from the sample/probe expression matrix and check that the optimum number of clusters is two (male and female), and generates a barplot of Y chromosome probe variance. A variance threshold can be selected, which reduces the number of Y chromosome probes used in the clustering step. 
The "massi.y" function extracts the Y chromosome probe values for each sample from the gene expression matrix and returns a matrix of sample/probe values that can be used directly with the "massi.cluster" function. The "massi.y" function also returns a barplot of Y chromosome probe variance showing the probe variance threshold as a .pdf file in the working directory.
The "massi.cluster" function clusters samples into two clustes using Y chromosome probe values. K-medoids clustering is performed using the partitioning around mediods (pam) method implimented in the "fpc" package.
The "massi.plot" function is for generating plots using the output from the massi.check and massi.cluster outputs.
The "massi.espresso" function combines all of the key features of MASSI and and performs the analysis using default values.

}
\author{
Sam Buckberry

Maintainer: Sam Buckberry <sam.buckberry@adelaide.edu.au> 
}
\references{
Christian Hennig (2013). fpc: Flexible procedures for clustering. R package version 2.1-5.
  http://CRAN.R-project.org/package=fpc
  
Gregory R. Warnes. Includes R source code and/or documentation contributed by: Ben Bolker, Lodewijk
  Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson,
  Steffen Moeller, Marc Schwartz and Bill Venables (2013). gplots: Various R programming tools for plotting
  data. R package version 2.11.0.1. http://CRAN.R-project.org/package=gplots
}

\keyword{ package }
\seealso{
\code{\link{massi.plot}, \link{massi.cluster}, \link{massi.espresso}, \link{massi.y}, \link{massi.check}}
}


