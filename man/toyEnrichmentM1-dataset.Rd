\name{toyEnrichmentM1-dataset}
\docType{data}
\alias{toyEnrichmentM1}
\alias{toyEnrichmentM1-dataset}

\title{Toy dataset with SNP data}

\description{
This data set gives an \code{\linkS4class{Enrichment}} object after,
\code{\link{initFiles}} and \code{\link{readEnrichment}} is ran.
The LD extend method used is "LD" (in \code{\link{readEnrichment}}: \code{extendMethod="ld"}).
Compute LD for all SNPs two by two.
Genome Build 37.3 (hg19).
}

\usage{
data(toyEnrichmentM1)
toyEnrichmentM1
}

\format{
See class \code{\linkS4class{Enrichment}} for details about the format.
}

% \source{}

% \references{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{datasets}
\keyword{toyEnrichmentM1}
\keyword{snpEnrichment}