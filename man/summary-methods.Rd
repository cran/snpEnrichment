\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,Chromosome-method}
% \alias{summary,Chromosome,ANY,ANY-method}
\alias{summary,Enrichment-method}
% \alias{summary,Enrichment,ANY,ANY-method}

\title{
Summary method (S4)
}

\description{
\code{\link{summary}} is a generic function used to produce result summaries of the results of various model fitting functions.
The function invokes particular \code{methods} which depend on the \code{class} of the first argument.
}

\usage{
\S4method{summary}{Enrichment}(object, extended = TRUE, complete = TRUE) 

\S4method{summary}{Chromosome}(object, extended = TRUE) 
}

\arguments{
    \item{object}{[Enrichment or Chromosome]: an object of class \code{\linkS4class{Enrichment}} or \code{\linkS4class{Chromosome}}.}
    \item{extended}{[logical]: \code{extended=TRUE} (default) to add \code{xSNP}.}
    \item{complete}{[logical]: \code{complete=TRUE} (default) to add chromosomes results.}
}

% \details{}

\value{
Return a \code{vector} (\code{extended=FALSE}) or a list of two \code{vector} (\code{extended=TRUE}) for classes \code{\linkS4class{Enrichment}} and  \code{\linkS4class{Chromosome}} (\code{complete=TRUE}).\cr
If \code{complete=TRUE}, return a \code{data.frame} (\code{extended=FALSE}) or a list of two \code{data.frame} (\code{extended=TRUE}) for class \code{\linkS4class{Enrichment}}.
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{snpEnrichment}
\keyword{Enrichment}
\keyword{summary}
\keyword{methods}

\examples{
data(toyM1)
summary(toyM1)
}
