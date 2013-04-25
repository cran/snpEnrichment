\name{writeLD}
\alias{writeLD}

\title{
Linkage Disequilibrium (LD) computation with PLINK
}

\description{
\code{\link{writeLD}} write a '.ld' file for each chromosomes which contains the LD (r^2).
}

\usage{
writeLD(pattern = "Chrom", snpInfoDir, signalFile, 
        ldThresh = 0.8, mc.cores = detectCores())
}

\arguments{
    \item{pattern}{[character]: character string containing a expression to be matched with all chromosomes files (e.g."Chrom" for  files which start by "Chrom" followed by the chromosome number).}
    \item{snpInfoDir}{[character]: character string naming a directory containing the reference data in a PLINK format (*.bed, *.bim and *.fam).}
    \item{signalFile}{[character]: the name of the signal file which the data are to be read from (2 columns: "SNP" and "PVALUE"). 
    Each row of the table appears as one line of the file.
    If it does not contain an \code{_absolute_} path, the file name is \code{_relative_} to the current working directory, \code{getwd}.
    The fields separator character have to be a space \code{" "} or a tabulation \code{"\\t"}.}
    \item{ldThresh}{[numeric]: threshold  value for LD calculation (PLINK).}
    \item{mc.cores}{[numeric]: The number of cores to use, i.e. at most how many child processes will be run simultaneously. 
    The option is initialized with \code{detectCores}. 
    Must be at least one, and parallelization requires at least two cores.}
}

% \details{}

\value{
One ".ld" file by chromosome is returned by \code{\link{writeLD}} in \code{snpInfoDir} directory.
}

\note{
The LD computation can take a long time depending on number of SNPs in \code{signalFile}.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{Enrichment}
\keyword{initFiles}
\keyword{LD}
\keyword{ld}
\keyword{linkage}
\keyword{Disequilibrium}
\keyword{PLINK}

\examples{
## Not run:
# signalFile <- "./extData/Signal/toySignal.txt"
# snpInfoDir <- "./extData/snpInfo/"
# writeLD(pattern = "Chrom", snpInfoDir, signalFile, 
#         ldThresh = 0.8, mc.cores = detectCores())
## End (Not run)
}
