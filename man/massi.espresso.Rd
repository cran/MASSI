\name{massi.espresso}
\alias{massi.espresso}

\title{massi.espresso}
\description{
This function combines all of the key features of MASSI and applies default values. 
}
\usage{
massi.espresso(exprs, y.probes, threshold = 3)
}

\arguments{
  \item{exprs}{
A sample/probe data.frame of normalized values with sample names as column names and probe id's as row names.
}
  \item{y.probes}{
A data.frame of probe id's in one column that match Y chromosome genes.
}
  \item{threshold}{
This specifies the threshold of probes to be included in the analysis. The threshold is determined by quantiles of probe variance (CV): 1=All probes, 2=Upper 75\%, 3=Upper 50\%, 4=Upper 25\%. It is highly recommended that probe CV plot generated using the massi.y function be inspected to inform threshold choice. Default=3.
}
}
\details{
See documentation for other other MASSI functions:
\code{\link{massi.plot}, \link{massi.cluster}, \link{massi.y}, \link{massi.check}}
}
\value{
k.medoids: the results of the k-medoids clustering of y chromosome probes.
"massi.results.txt" is a tab-delimited file saved to the current working directory.
"MASSI.figures.pdf" is a collection of figures generated using the massi.plot function and saved to the current working directory.
}
\references{
Christian Hennig (2013). fpc: Flexible procedures for clustering. R package version 2.1-5.
  http://CRAN.R-project.org/package=fpc
  
Gregory R. Warnes. Includes R source code and/or documentation contributed by: Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley,
Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2013).

gplots: Various R programming tools for plotting data. R package version 2.11.0.1.
http://CRAN.R-project.org/package=gplots

}
\author{
Sam Buckberry
}

\seealso{
\code{\link{massi.plot}, \link{massi.cluster}, \link{massi.y}, \link{massi.check}}

}
\examples{
\dontrun{
data(MASSI.test.dataset, MASSI.test.probes)
massi.espresso(MASSI.test.dataset, MASSI.test.probes)
}
}

