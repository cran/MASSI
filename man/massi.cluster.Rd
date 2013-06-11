\name{massi.cluster}
\alias{massi.cluster}

\title{massi.cluster}
\description{
Clusters samples into two clustes using Y chromosome probe values. K-medoids clustering is performed using the partitioning around mediods (pam) method implimented in the "fpc" package. The cluster with the highest probes values is determined to be the cluster of "male" samples and the cluster with the lowest values as female samples. A tab-delimited file is output to the working directory which includes the sample id, mean probe value for each sample, the probe standard deviation for each sample and the sex assigned to each sample.
}
\usage{
massi.cluster(y.subset.values)
}

\arguments{
  \item{y.subset.values}{
This object is created using the "massi.check" function. These are probe values from the Y chromosome probes selected in the "massi.check" step. This data.frame has the sample names as column names and the row names as probe id's.  
}
}

\value{
k.medoids: the results of the k-medoids clustering of y chromosome probes.
"massi.results.txt" is a tab-delimited file saved to the current working directory.
}
\references{
  Christian Hennig (2013). fpc: Flexible procedures for clustering. R package version 2.1-5.
  http://CRAN.R-project.org/package=fpc
}
\author{
Sam Buckberry
}

\seealso{
\code{\link{massi.plot}, \link{massi.y}, \link{massi.espresso}, \link{massi.check}}
}

\examples{
\dontrun{
data(MASSI.test.dataset, MASSI.test.probes)
massi.cluster(y.subset.values)
}
}

