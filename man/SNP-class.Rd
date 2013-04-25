\name{SNP-class}
\docType{class}
\alias{SNP-class}
\alias{SNP}
\alias{[,SNP-method}
\alias{[<-,SNP-method}
\alias{show,SNP-method}
% \alias{newSNP}
% \alias{newSNP-methods}
% \alias{newSNP,ANY,ANY,ANY,ANY,ANY,ANY-method}

\title{Class \code{"\linkS4class{SNP}"}}

\description{
This class is defined to summarize the enrichment analysis. 
It's a part of \code{\linkS4class{Chromosome}} and \code{\linkS4class{Enrichment}} classes.
}

% \section{Objects from the Class}{}

\section{Slots}{
    \describe{
        \item{List}{[vector(character)]: a list of SNPs used to compute enrichment (e.g. eSNP or xSNP).}
        \item{Table}{[matrix]: Contingency table with SNPs (columns) and P-Values from signal (rows).}
        \item{EnrichmentRatio}{[numeric]: Enrichment Ratio is computed on the contingency table (\code{Table} slot).}
        \item{Z}{[numeric]: A statistic computed from \code{EnrichmentRatio} and resampling results.}
        \item{PValue}{[numeric]: P-Value associated with the statistic \code{Z}.}
        \item{Resampling}{[matrix]: A matrix with by row, the contingency table and the odds ratio for each resampling.}
    }
}

\section{Methods}{
    \describe{
        % \item{newSNP(List, Table, EnrichmentRatio, Z, PValue, Resampling): }{Generate and initialize a new \code{\linkS4class{SNP}} object.}
        \item{object["slotName"]: }{Get the value of the field \code{slotName}.}
        \item{object["slotName"]<-value: }{Set \code{value} to the field \code{slotName}.}
        \item{show(object): }{Return the formatted values of \code{\linkS4class{SNP}} object.}
    }
}

\note{
\code{\linkS4class{SNP}} object is not intended to be use directly by user. 
It is a part of the \code{\linkS4class{Enrichment}} and \code{\linkS4class{Chromosome}}object.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{classes}
\keyword{class}
\keyword{snp}

% \examples{}

