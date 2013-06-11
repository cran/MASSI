\name{y.probes}
\alias{y.probes}
\docType{data}
\title{
Y chromosome probes
}
\description{
lists of y chromosome probes corresponding to commercial microarray platforms.
}
\usage{data(y.probes)}
\format{
  A data frame with 171 observations on the following 8 variables.
  \describe{
    \item{\code{Illumina.Human.HT.12.probe}}{a character vector}
    \item{\code{Illumina.HumanWG.6.v1.probe}}{a character vector}
    \item{\code{Illumina.HumanWG.6.v2.probe}}{a character vector}
    \item{\code{Illumina.HumanWG.6.v3.probe}}{a character vector}
    \item{\code{Affy.HG.U133.PLUS.2.probeset}}{a character vector}
    \item{\code{Affy.HuGene.1_0.st.v1.probeset}}{a character vector}
    \item{\code{Affy.MoGene.probeset}}{a character vector}
    \item{\code{Affy.MoEx.probeset}}{a character vector}
  }
}
\details{
Microarray probes were mapped to human or mouse genomes to establish which probes represent Y chromsome gene using the methods detailed at http://jan2013.archive.ensembl.org/info/docs/microarray_probe_set_mapping.html.
}

\examples{
data(y.probes)
## maybe str(y.probes) ; plot(y.probes) ...
}
\keyword{datasets}
