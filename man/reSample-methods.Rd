\name{reSample}
\alias{reSample}
\docType{methods}
\alias{reSample-methods}
\alias{reSample,Chromosome-method}
% \alias{reSample,Chromosome,ANY,ANY,ANY,ANY,ANY-method}
% \alias{reSample,Chromosome,ANY-method}
\alias{reSample,Enrichment-method}
% \alias{reSample,Enrichment,ANY,ANY,ANY,ANY,ANY-method}
% \alias{reSample,Enrichment,ANY-method}
\alias{reSample,ANY-method}

\title{
Compute enrichment analysis on an \code{\linkS4class{Enrichment}} object
}

\description{
After \code{\link{initFiles}} and \code{\link{readEnrichment}} has been run. 
\code{\link{reSample}} computes a statistic value and a p-value for each chromosomes and for the whole genome.
}

\usage{
reSample(object, nSample = 100, sigThresh = 0.05, 
         MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
         extendMethod = "ld", mc.cores = detectCores())
}

\arguments{
    \item{object}{[Enrichment or Chromosome]: an object to be updated. It is intended, an object returned by the \code{\link{readEnrichment}} function.}
    \item{nSample}{[numeric]: the number of resampling done by \code{\link{reSample}} for p-values computation (minimum is 100).}
    \item{sigThresh}{[numeric]: statistical threshold for signal (e.g. \code{sigThresh = 0.05} for a given GWAS signal) used to compute an Enrichment Ratio.}
    \item{MAFpool}{[vector(numeric)]: either a numeric vector giving the breaks points of intervals into which SNP's MAF (Minor Allele Frequency) is to be split.}
    \item{extendMethod}{[character]: character string naming the method to be used to extend the list of eSNPs. 
    Default is "ld", the list of eSNPs is extended with all SNPs in LD (r^2 threshold, defined in \code{\link{initFiles}}).
    Otherwise, "block", extends the list of eSNPs by defining LD blocks with LD measures and positions of SNPs on chromosomes.}
    \item{mc.cores}{[numeric]: The number of cores to use, i.e. at most how many child processes will be run simultaneously. 
    The option is initialized with \code{detectCores}. 
    Must be at least one, and parallelization requires at least two cores.}
}

% \details{}

\value{
Return the object given in argument, updated by the resampling results.
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{snpEnrichment}
\keyword{Enrichment}
\keyword{reSample}
\keyword{methods}

\examples{
# data(toyM1)
# reSample(object = toyM1, 
         # nSample = 10, 
         # sigThresh = 0.05, 
         # MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), 
         # extendMethod = "ld", mc.cores = detectCores())
# toyM1
}
