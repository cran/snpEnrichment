\name{excludeSNP}
\alias{excludeSNP}
\docType{methods}
\alias{excludeSNP-methods}
\alias{excludeSNP,Enrichment-method}
\alias{excludeSNP,ANY-method}
% \alias{excludeSNP,Enrichment,ANY-method}
% \alias{excludeSNP,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method}

\title{
Exclude SNPs from Enrichment analysis
}

\description{
Remove a specify set of SNPs and compute a new enrichment analysis.
}

\usage{
excludeSNP(object, excludeFile, nSample = 100, sigThresh = 0.05, 
           MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), 
           extendMethod = "ld", mc.cores = detectCores())
}

\arguments{
    \item{object}{[Enrichment]: an \code{\linkS4class{Enrichment}} object filled by \code{\link{reSample}}.}
    \item{excludeFile}{[vector(character)]: a list of SNPs to remove from a previous enrichment analysis. A path to a file which the first column are the SNPs.}
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
Return the object given in argument where lists of SNPs are updated by removing SNPs in \code{excludeFile}.
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

\keyword{enrichment}
\keyword{excludeSNP}
\keyword{methods}

\examples{
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
}
