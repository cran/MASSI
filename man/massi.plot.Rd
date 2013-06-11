\name{massi.plot}
\alias{massi.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{massi.plot}
\description{
A function for generating plots from the massi.check and massi.cluster outputs
}
\usage{
massi.plot(y.values = y.values, y.subset.values = y.subset.values,
  massi.output = massi.output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.values}{
This object is created using the "massi.check" function. These are all of the Y chromosome probes values for each sample.
}
  \item{y.subset.values}{
This object is created using the "massi.check" function. These are the subset of y.values object output from the "massi.check" step. This data.frame has the sample names as column names and the row names as probe id's.
}
  \item{massi.output}{
A matrix of values generated using the massi.cluster function. 
}
}
\details{
This function returns a pdf file "MASSI.figures.pdf" in the current working directory. 
}

\references{
Gregory R. Warnes. Includes R source code and/or documentation contributed by: Ben Bolker, Lodewijk
  Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson,
  Steffen Moeller, Marc Schwartz and Bill Venables (2013). gplots: Various R programming tools for plotting
  data. R package version 2.11.0.1. http://CRAN.R-project.org/package=gplots
}
\author{
Sam Buckberry
}

\seealso{
\code{\link{massi.cluster}, \link{massi.y}, \link{massi.espresso}, \link{massi.check}}
}
\examples{
\dontrun{
massi.plot(y.values=y.values, y.subset.values=y.subset.values,
                       massi.output=massi.output)
}
}

