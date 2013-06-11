\name{massi.check}
\alias{massi.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{massi.check}
%%  ~~function to do ... ~~

\description{This function determines if the optimum number of clusters (k) in the dataset is 2; which is highly desirable for downstream functions in MASSI. If the optimum number of clusters is 2, then the message "Good, Optimum number of clusters is 2..." will be printed. If the optimum number of clusters is >2, then the message "Bad, Optimum number of clusters is >2, try reducing the Y chromosome probe list by setting a higher threshold" will be printed. 
}

\usage{
massi.check(exprs, y.probes, threshold = 3)
}

\arguments{
  \item{exprs}{
A sample/probe data.frame of normalized values with sample names as column names and probe id's as row names. 
}
  \item{y.probes}{
A data.frame of probe id's in one column that match Y chromosome genes
}
  \item{threshold}{
This specifies the threshold of probes to be included in the analysis. The threshold is determined by quantiles of probe variance (CV): 1=All probes, 2=Upper 75\%, 3=Upper 50\%, 4=Upper 25\%. It is highly recommended that probe CV plot generated using the massi.y function be inspected to inform threshold choice. Default=3. If the total number of Y chromosome probes in the dataset is less than 10, consider using a threshold of 1 (all probes).
}
}
\details{
The optimum number of clusters in the dataset is estimated using the average silhouette width of clusters as implimented in the 'fpc' package.
The returned matrix "y.subset.values" is used by the "massi.cluster" function.
}
\value{
y.values: a matrix of Y chromosome probe values
y.subset.values: a subset matrix of y.values for probes above selected threshold.
}
\references{
  Christian Hennig (2013). fpc: Flexible procedures for clustering. R package version 2.1-5.
  http://CRAN.R-project.org/package=fpc
}
\author{
Sam Buckberry
}

\seealso{
\code{\link{massi.plot}, \link{massi.cluster}, \link{massi.espresso}, \link{massi.y}}
}

\examples{
\dontrun{
data(MASSI.test.dataset, MASSI.test.probes)
massi.check(exprs=MASSI.test.dataset, y.probes=MASSI.test.probes)
}
}


