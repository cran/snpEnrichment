\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,Enrichment-method}
\alias{plot,Enrichment,ANY-method}

\title{
Plot method (S4) for \code{\linkS4class{Enrichment}} object
}

\description{
\code{\link{plot}} is a generic function for plotting of R objects.
The function invokes particular \code{methods} which depend on the \code{class} of the first argument.
}

\usage{
\S4method{plot}{Enrichment}(x, onlyGenome = FALSE, chr = NULL, 
           types = c("eSNP", "xSNP"), ...)
}

\arguments{
    \item{x}{[Enrichment]: an object of class \code{\linkS4class{Enrichment}} which the Z statistics have to be drawn.}
    \item{onlyGenome}{[character]: if \code{TRUE} (default \code{onlyGenome=FALSE}) plot Z statistics for genome only.}
    \item{chr}{[vector(numeric)]: \code{chr=NULL} (default) plot Z statistics for the whole genome and all chromosomes or the selected chromosomes.}
    \item{types}{[vector(character)]: plot the selected analysis for "eSNP" and/or "xSNP".}
    \item{...}{[any]: Arguments to be passed to methods, such as graphical parameters (see \code{par})}
}

% \details{}

% \value{}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{EnrichSNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{plot}
\keyword{methods}

\examples{
\dontrun{
data(toyEnrichmentM1)
reSample(toyEnrichmentM1, 10)
plot(toyEnrichmentM1)}
}
