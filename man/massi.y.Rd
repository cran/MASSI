\name{massi.y}
\alias{massi.y}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{massi.y}
\description{
This function extracts the Y chromosome probe values for each sample from the gene expression matrix.
}
\usage{
massi.y(exprs, y.probes)
}

\arguments{
  \item{exprs}{
A gene expression matrix in the data.frame format.
Samples names must be column names and probe id's as row names
}
  \item{y.probes}{
A list of Y chromosome probe id's in the data.frame format.
This matrix must consist of only one column.
}
}
\details{
This function extracts expression values for probes that correspond to Y
chromosome genes from the sample/probe expression matrix.
At this step, you can use your own list of Y chromosome probes
or using the probe lists included with the MASSI package that correspond
to commonly used commercial microarray platforms.
This function produces a plot the Y chromosome probe variation across samples and
indicates quantiles used as thresholds and is saved as a pdf file in the
current working directory.
This plot should be inspected prior to running the massi.check function.

}
\value{
massi.yplot.pdf
}

\author{
Sam Buckberry
}

\seealso{
\code{\link{massi.plot}, \link{massi.cluster}, \link{massi.espresso}, \link{massi.check}}
}

\examples{
\dontrun{
data(MASSI.test.dataset, MASSI.test.probes)
massi.y(exprs=MASSI.test.dataset, y.probes=MASSI.test.probes)
}
}

