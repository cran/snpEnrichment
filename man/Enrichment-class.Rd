\name{Enrichment-class}
\docType{class}
\alias{Enrichment-class}
\alias{Enrichment}
\alias{[,Enrichment-method}
\alias{[<-,Enrichment-method}
% \alias{computeER,Enrichment-method}
% \alias{excludeSNP,Enrichment-method}
% \alias{reSample,Enrichment-method}
% \alias{reset,Enrichment-method}
\alias{show,Enrichment-method}
% \alias{summary,Enrichment-method}
\alias{enrichment}
\alias{enrichment-methods}
% \alias{enrichment,ANY,ANY,ANY,ANY,ANY-method}
\alias{enrichment,ANY-method}

\title{Class \code{\linkS4class{Enrichment}}}

\description{
This class is defined to summarize the enrichment analysis on each chromosomes and the whole genome.
}

\section{Objects from the Class}{
    \code{\link{enrichment}} is defined to build an object of class \code{\linkS4class{Enrichment}} in order to compute an enrichment analysis.
    \code{Enrichment} is the object containing the results for all \code{\linkS4class{Chromosome}} object and for the whole genome.\cr

    When an \code{\linkS4class{Enrichment}} object is created, it contains the original signal (e.g. GWAS), a list of SNPs. 
    All the others slots are "empty". After \code{\link{reSample}} is ran on an \code{\linkS4class{Enrichment}} object, 
    the slots: Table, EnrichmentRatio, Z, PValue and Resampling are filled.

    Note that if \code{\link{reSample}} is executed on an \code{\linkS4class{Enrichment}} every new resampling
    is added to the original ones, pre-existing statistics are erased and computed again with the new resampling set.
}

\section{Slots}{
    \describe{
        \item{Signal}{[data.frame]: a three columns \code{data.frame}: "SNP", "PVALUE" and "IN" (e.g. GWAS). "IN" is computed during the reading step and gives informations about which SNPs are kept for the enrichment analysis.}
        \item{Parameters}{[list]: each parameters used for the reading or resampling step are stored in this slot.}
        \item{eSNP, xSNP}{[SNP]: contain a \code{\linkS4class{SNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for a list of SNPs (\code{eSNP}) and an extended one (\code{xSNP}).}
        % \item{xSNP}{[SNP]: contains a \code{\linkS4class{SNP}} object (whith slots: SNP, Table, EnrichmentRatio, Z, PValue and Resampling) for an extended list of SNPs.}
        \item{Chromosomes}{[list(Chromosome)]: a list of 22 \code{\linkS4class{Chromosome}} objects.}
    }
}

\section{Extends}{
Class \code{\linkS4class{Chromosome}}, directly.\cr
Class \code{\linkS4class{SNP}}, directly.\cr
}

\section{Methods}{
    \describe{
        \item{enrichment(Signal, Parameters, eSNP, xSNP, Chromosomes):}{Generate and initialize a new \code{\linkS4class{Enrichment}} object.}
        \item{object["slotName"]:}{Get the value of the field \code{slotName}.}
        \item{object["slotName"]<-value:}{Set \code{value} to the field \code{slotName}.}
        \item{show(object):}{Return the formatted values of \code{\linkS4class{Enrichment}} object.}
        % \item{computeER(object, sigThresh):}{Compute Enrichment Ratio for \code{\linkS4class{Enrichment}} object (internal use only).}
        \item{reSample(object, nSample, sigThresh, MAFpool, extendMethod, mc.cores):}{Compute P-Values based upon a resampling of SNPs (\code{eSNP} and \code{xSNP}) and update an \code{\linkS4class{Enrichment}} object.}
        \item{excludeSNP(object, excludeFile, nSample, sigThresh, MAFpool, extendMethod, mc.cores):}{Excludes SNPs given in \code{excludeFile} from the original list of eSNPs (xSNPs). Then a new enrichment analysis is computed.}
        \item{summary(object, extended, complete):}{Return a list of \code{data.frame} which summarize enrichment analysis for \code{eSNP} and/or \code{xSNP}.
            \code{extended=TRUE} (default) to add \code{xSNP}.
            \code{complete=TRUE} (default) to add chromosomes results.}
        % \item{reset(object, "slotName"):}{Reset the field \code{slotName}.}
    }
}

% \note{}

\seealso{
Overview : \code{\link{snpEnrichment-package}} \cr
Classes : \code{\linkS4class{Enrichment}}, \code{\linkS4class{Chromosome}}, \code{\linkS4class{SNP}} \cr
Methods : \code{\link{reset}}, \code{\link{reSample}}, \code{\link{excludeSNP}}, \code{\link{compareEnrichment}}, \cr
            \code{\link{enrichment}}, \code{\link{is.enrichment}}, \code{\link{chromosome}}, \code{\link{is.chromosome}} \cr
Functions : \code{\link{initFiles}}, \code{\link{writeLD}}, \code{\link{readEnrichment}}
}

\keyword{classes}
\keyword{class}
\keyword{enrichment}
\keyword{chromosome}
\keyword{snp}

\examples{
data(toyM1)
toyEnrich <- enrichment()
show(toyEnrich)

toyEnrich["Signal"] <- toyM1["Signal"]
head(toyEnrich["Signal"])

toyEnrich <- enrichment(Signal = toyM1["Signal"], eSNP = toyM1["eSNP"])
toyEnrich <- enrichment(Signal = toyM1["Signal"])

# reSample(object = toyM1, 
#          nSample = 10, 
#          sigThresh = 0.05, 
#          MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5),
#          extendMethod = "ld",
#          mc.cores = detectCores())
summary(toyM1)

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
# summary(toyM1_exclude)
}
