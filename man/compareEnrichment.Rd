\name{compareEnrichment}
\alias{compareEnrichment}
% \docType{methods}
\alias{compareEnrichment-methods}
% \alias{compareEnrichment,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}
\alias{compareEnrichment,ANY-method}

\title{
Compare enrichment analysis between two SNPs list
}

\description{
Compare the enrichment analysis between two set of SNPs. 
\code{\link{compareEnrichment}} compare two files with list of SNPs or two \code{\linkS4class{Enrichment}} objects.
}

\usage{
compareEnrichment(list1, list2, pattern = "Chrom", signalFile, 
                  transcriptFile = FALSE, snpInfoDir, 
                  distThresh = 1000, LD = FALSE, 
                  ldThresh = 0.8, nSample = 100, 
                  sigThresh = 0.05, 
                  MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
                  extendMethod = "ld",
                  mc.cores = detectCores())
}

\arguments{
    \item{list1, list2}{[character or Enrichment]: character string naming a directory containing a list of SNPs for each chromosomes or an \code{\linkS4class{Enrichment}} object (fully filled).}
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{signalFile}{[character]: the name of the signal file which the data are to be read from (2 columns: "SNP" and "PVALUE"). 
    Each row of the table appears as one line of the file.
    If it does not contain an \code{_absolute_} path, the file name is \code{_relative_} to the current working directory, \code{getwd}.
    The fields separator character have to be a space \code{" "} or a tabulation \code{"\\t"}.}
    \item{transcriptFile}{[character or data.frame]: character string naming a file or a \code{data.frame} with four columns: Chromomosome, trancript's name, Starting and Ending positions. 
    \code{data(trancript)} can be use as parameters. Default is \code{FALSE}.}
    \item{snpInfoDir}{[character]: character string naming a directory containing the reference data in a PLINK format (*.bed, *.bim and *.fam).}
    \item{distThresh}{[numeric]: maximal distance (kb) between SNP and gene. \code{distThresh} is used if \code{transcriptFile} is set.}
    \item{LD}{[logical]: \code{LD=TRUE} (default is \code{FALSE}) compute LD with PLINK.
    Note that, this setting can increase the computation's time, depending on number of SNPs in the signal file.}
    \item{ldThresh}{[numeric]: threshold  value for LD calculation (PLINK).}
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
Return a \code{list} of two \code{data.frame} which summarize for \code{eSNP} and/or \code{xSNP} the enrichment analyses and the comparison between them.
}

% \note{
% Still in development.
% }

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{plot}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{enrichment}
\keyword{chromosome}
\keyword{compareEnrichment}

% \examples{}
