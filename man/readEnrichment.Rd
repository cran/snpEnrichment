\name{readEnrichment}
\alias{readEnrichment}

\title{
Read and create EnrichmentRatio object
}

\description{
Read files created by \code{\link{initFiles}} and create an \code{\linkS4class{Enrichment}} object.
}

\usage{
readEnrichment(pattern = "Chrom", signalFile, 
               transcriptFile = "FALSE", snpListDir, 
               snpInfoDir, distThresh = 1000, 
               sigThresh = 0.05, LD = "FALSE", 
               extendMethod = "ld", mc.cores = detectCores())
}

\arguments{
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{signalFile}{[character]: the name of the signal file which the data are to be read from (2 columns: "SNP" and "PVALUE"). 
    Each row of the table appears as one line of the file.
    If it does not contain an \code{_absolute_} path, the file name is \code{_relative_} to the current working directory, \code{getwd}.
    The fields separator character have to be a space \code{" "} or a tabulation \code{"\\t"}.}
    \item{transcriptFile}{[character or data.frame]: character string naming a file or a \code{data.frame} with four columns: Chromomosome, trancript's name, Starting and Ending positions. 
    \code{data(trancript)} can be use as parameters. Default is \code{FALSE}.}
    \item{snpListDir}{[character]: character string naming a directory containing a list of SNPs for each chromosomes.}
    \item{snpInfoDir}{[character]: character string naming a directory containing the reference data in a PLINK format (*.bed, *.bim and *.fam).}
    \item{distThresh}{[numeric]: maximal distance (kb) between SNP and gene. \code{distThresh} is used if \code{transcriptFile} is set.}
    \item{sigThresh}{[numeric]: statistical threshold for signal (e.g. \code{sigThresh = 0.05} for a given GWAS signal) used to compute an Enrichment Ratio.}
    \item{LD}{[logical]: \code{LD=TRUE} (default is \code{FALSE}) compute LD with PLINK.
    Note that, this setting can increase the computation's time, depending on number of SNPs in the signal file.}
    \item{extendMethod}{[character]: character string naming the method to be used to extend the list of eSNPs. 
    Default is "ld", the list of eSNPs is extended with all SNPs in LD (r^2 threshold, defined in \code{\link{initFiles}}).
    Otherwise, "block", extends the list of eSNPs by defining LD blocks with LD measures and positions of SNPs on chromosomes.}
    \item{mc.cores}{[numeric]: The number of cores to use, i.e. at most how many child processes will be run simultaneously. 
    The option is initialized with \code{detectCores}. 
    Must be at least one, and parallelization requires at least two cores.}
}

% \details{}

\value{
Return an object of class \code{\linkS4class{Enrichment}} partly filled.
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{enrichment}
\keyword{initFiles}
\keyword{writeLD}
\keyword{initialize}
\keyword{readEnrichment}

\examples{
## Not run:
# snpListDir <- "./extData/List/"
# signalFile <- "./extData/Signal/toySignal.txt"
# excludeFile <- "./extData/Exclude/toyExclude.txt"
# snpInfoDir <- "./extData/snpInfo/"
# data(transcript)
# transcriptFile <- transcript

# toy_M1 <- readEnrichment(pattern = "Chrom", signalFile, 
#                          transcriptFile, snpListDir, 
#                          snpInfoDir, distThresh = 1000, 
#                          sigThresh = 0.05, LD = "FALSE", 
#                          extendMethod = "ld", 
#                          mc.cores = detectCores())
# toy_M1
## End (Not run)
}
