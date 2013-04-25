\name{snpEnrichment-package}
\alias{snpEnrichment-package}
\docType{package}
\alias{snpEnrichment}

\title{~ Overview: SNPs enrichment analysis ~}

\description{
Implements classes and methods for large-scale SNP enrichment analysis (e.g. SNPs associated with genes expression and a GWAS signal).
}

\details{\tabular{ll}{
Package: \tab snpEnrichment\cr
Title: \tab SNPs enrichment analysis\cr
Author: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
Contributor: \tab Loic Yengo <loic.yengo at good.ibl.fr>\cr
Maintainer: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
% Description: \tab Implements classes and methods for large-scale SNP enrichment analysis \cr (e.g. SNPs associated with genes expression and a GWAS signal) \cr
License: \tab GPL (>= 2)\cr
Depends: \tab R (>= 2.15.0), methods, parallel\cr
URL: \tab http://www-good.ibl.fr/\cr
Encoding: \tab UTF-8\cr
}}

\note{
In order to use all options, PLINK must be installed. PLINK is used to compute the Linkage Disequilibrium and get the Minor Allele Frequency from a reference dataset like 1Kg Genome.
}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods :  \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{package}
\keyword{snpEnrichment}
\keyword{Enrichment}

\examples{
#######################
### 1. Data Preparation
## Not run:
# snpInfoDir <- "./extdata/snpInfo/"
# signalFile <- "./extdata/Signal/toySignal.txt"
# initFiles(pattern = "Chrom", snpInfoDir, signalFile, 
#           ldThresh = 0.8, LD = TRUE, mc.cores = detectCores())
## End (Not run)

## OR
## Not run:
# snpInfoDir <- "./extdata/snpInfo/"
# signalFile <- "./extdata/Signal/toySignal.txt"
# initFiles(pattern = "Chrom", snpInfoDir, signalFile, 
#           ldThresh = 0.8, LD = FALSE, mc.cores = detectCores())
# writeLD(pattern = "Chrom", snpInfoDir, signalFile, 
#         ldThresh = 0.8, mc.cores = detectCores())
## End (Not run)


###################
### 2. Reading data
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


########################
### 3. Computing results
# data(toyM1)
# reSample(object = toyM1, 
         # nSample = 10, 
         # sigThresh = 0.05, 
         # MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), 
         # extendMethod = "ld", 
         # mc.cores = detectCores())
# toyM1


#######################
### 5. Further analysis
## Exclude SNP from original list.
# data(toyM1) # data(toyM2)
# excludeFile <- c(
    # "rs1561385", "rs7792796", "rs2514670", "rs9641913", "rs8184976",
    # "rs17582442", "rs7690663", "rs4940941", "rs7069561", "rs540218",
    # "rs17315714", "rs17795475", "rs7171423", "rs2392927", "rs12593911",
    # "rs4150477", "rs11608342", "rs16998578", "rs4299828", "rs915865",
    # "rs10976361", "rs7863276", "rs16908503", "rs544845", "rs1473462",
    # "rs4757541", "rs7640480", "rs7121036", "rs6803546", "rs10851981",
    # "rs4724502", "rs9540053", "rs10935849", "rs11193005", "rs6566417",
    # "rs1693294", "rs12759271", "rs17718970", "rs4774717", "rs455839",
    # "rs942278", "rs6545708", "rs7557832", "rs1498356", "rs11083318",
    # "rs9595937", "rs1561476", "rs12188654", "rs2048839", "rs4689801"
# )

# toyM1_exclude <- excludeSNP(toyM1, 
                            # excludeFile, 
                            # nSample = 10, 
                            # sigThresh = 0.05, 
                            # MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
                            # extendMethod = "ld", 
                            # mc.cores = detectCores())
# toyM1_exclude


####################
### 4. Watch results
# show(toyM1)
# summary(toyM1)

# show(toyM1_exclude)
# summary(toyM1_exclude)
}
