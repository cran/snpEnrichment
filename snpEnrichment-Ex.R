pkgname <- "snpEnrichment"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('snpEnrichment')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Chromosome-class")
### * Chromosome-class

flush(stderr()); flush(stdout())

### Name: Chromosome-class
### Title: Class 'Chromosome'
### Aliases: Chromosome-class Chromosome [,Chromosome-method
###   [<-,Chromosome-method show,Chromosome-method chromosome
###   chromosome-methods chromosome,ANY-method
### Keywords: classes class methods method chromosome snp

### ** Examples

Data <- structure(
    list(
        SNP = c("rs4970420", "rs3766178", 
                "rs10910030", "rs10910036", 
                "rs2234167", "rs6661861"), 
        PVALUE = c(0.9244, 0.167, 0.01177, 0.4267, 0.9728, 0.4063), 
        CHR = c(1, 1, 1, 1, 1, 1),
        POS = c(1106473, 1478180, 2035684, 2183754, 2494330, 3043121), 
        MAF = c(0.2149, 0.3117, 0.374, 0.3753, 0.1565, 0.06101), 
        eSNP = c(0, 1, 1, 0, 0, 0), 
        xSNP = c(0, 1, 1, 0, 0, 0)
    ), 
    .Names = c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP", "xSNP"), 
    row.names = c("rs4970420", "rs3766178", 
                  "rs10910030", "rs10910036", 
                  "rs2234167", "rs6661861"),
class = "data.frame")

toyChr <- chromosome(Data = Data)
show(toyChr)
toyChr

toyChr <- chromosome()
toyChr["Data"] <- Data
toyChr

is.chromosome(toyChr) # TRUE



cleanEx()
nameEx("Enrichment-class")
### * Enrichment-class

flush(stderr()); flush(stdout())

### Name: Enrichment-class
### Title: Class 'Enrichment'
### Aliases: Enrichment-class Enrichment [,Enrichment-method
###   [<-,Enrichment-method show,Enrichment-method enrichment
###   enrichment-methods enrichment,ANY-method
### Keywords: classes class enrichment chromosome snp

### ** Examples

data(toyM1)
toyEnrich <- enrichment()
show(toyEnrich)

toyEnrich["Loss"] <- toyM1["Loss"]
toyEnrich["Loss"]

toyEnrich <- enrichment(Loss = toyM1["Loss"], eSNP = toyM1["eSNP"])
toyEnrich <- enrichment(Loss = toyM1["Loss"])

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



cleanEx()
nameEx("GC")
### * GC

flush(stderr()); flush(stdout())

### Name: GC
### Title: Full Garbage Collection
### Aliases: GC
### Keywords: GC garbage

### ** Examples

GC() #- do it now
x <- integer(100000); for(i in 1:18) x <- c(x,i)
GC(TRUE)
GC(reset=TRUE)



cleanEx()
nameEx("excludeSNP-methods")
### * excludeSNP-methods

flush(stderr()); flush(stdout())

### Name: excludeSNP
### Title: Exclude SNPs from Enrichment analysis
### Aliases: excludeSNP excludeSNP-methods excludeSNP,Enrichment-method
###   excludeSNP,ANY-method
### Keywords: enrichment excludeSNP methods

### ** Examples

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



cleanEx()
nameEx("initFiles")
### * initFiles

flush(stderr()); flush(stdout())

### Name: initFiles
### Title: Initialize files for enrichment analysis
### Aliases: initFiles
### Keywords: Enrichment initFiles writeLD initialize

### ** Examples

## Not run:
# snpInfoDir <- "./extdata/snpInfo/"
# signalFile <- "./extdata/Signal/toySignal.txt"
# snpListDir <- "./extData/List/"
# initFiles(pattern = "Chrom", snpInfoDir, snpListDir, signalFile, 
#           ldThresh = 0.8, LD = FALSE, mc.cores = detectCores())
## End (Not run)



cleanEx()
nameEx("is.chromosome")
### * is.chromosome

flush(stderr()); flush(stdout())

### Name: is.chromosome
### Title: Is an Chromosome object
### Aliases: is.chromosome is.chromosome-methods is.chromosome,ANY-method
### Keywords: chromosome is

### ** Examples

a <- chromosome()
c <- chromosome()
is.chromosome(list())                # FALSE
is.chromosome(1)                     # FALSE
is.chromosome(a)                     # TRUE
is.chromosome(c(a, c))               # TRUE TRUE
is.chromosome(list(a, b = "char"))   # TRUE FALSE
is.chromosome(c(a, b = list(12, c))) # TRUE FALSE TRUE



cleanEx()
nameEx("is.enrichment")
### * is.enrichment

flush(stderr()); flush(stdout())

### Name: is.enrichment
### Title: Is an Enrichment object
### Aliases: is.enrichment is.enrichment-methods is.enrichment,ANY-method
### Keywords: enrichment is

### ** Examples

a <- enrichment()
c <- enrichment()
is.enrichment(list())                # FALSE
is.enrichment(1)                     # FALSE
is.enrichment(a)                     # TRUE
is.enrichment(c(a, c))               # TRUE TRUE
is.enrichment(list(a, b = "char"))   # TRUE FALSE
is.enrichment(c(a, b = list(12, c))) # TRUE FALSE TRUE



cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

### Name: plot-methods
### Title: Plot method (S4) for 'Enrichment' object
### Aliases: plot plot-methods plot,Enrichment-method
###   plot,Enrichment,ANY-method
### Keywords: snpEnrichment Enrichment plot methods

### ** Examples

# data(toyM1)
# reSample(toyM1)
# plot(toyM1)



cleanEx()
nameEx("reSample-methods")
### * reSample-methods

flush(stderr()); flush(stdout())

### Name: reSample
### Title: Compute enrichment analysis on an 'Enrichment' object
### Aliases: reSample reSample-methods reSample,Chromosome-method
###   reSample,Enrichment-method reSample,ANY-method
### Keywords: snpEnrichment Enrichment reSample methods

### ** Examples

# data(toyM1)
# reSample(object = toyM1, 
         # nSample = 10, 
         # sigThresh = 0.05, 
         # MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), 
         # extendMethod = "ld", mc.cores = detectCores())
# toyM1



cleanEx()
nameEx("readEnrichment")
### * readEnrichment

flush(stderr()); flush(stdout())

### Name: readEnrichment
### Title: Read and create EnrichmentRatio object
### Aliases: readEnrichment
### Keywords: enrichment initFiles writeLD initialize readEnrichment

### ** Examples

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



cleanEx()
nameEx("snpEnrichment-package")
### * snpEnrichment-package

flush(stderr()); flush(stdout())

### Name: snpEnrichment-package
### Title: ~ Overview: SNPs enrichment analysis ~
### Aliases: snpEnrichment-package snpEnrichment
### Keywords: package snpEnrichment Enrichment

### ** Examples

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
#         ldThresh = 0.8, onlySignal = TRUE, 
#         mc.cores = detectCores())
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



cleanEx()
nameEx("summary-methods")
### * summary-methods

flush(stderr()); flush(stdout())

### Name: summary-methods
### Title: Summary method (S4)
### Aliases: summary summary-methods summary,Chromosome-method
###   summary,Enrichment-method
### Keywords: snpEnrichment Enrichment summary methods

### ** Examples

data(toyM1)
summary(toyM1)



cleanEx()
nameEx("writeLD")
### * writeLD

flush(stderr()); flush(stdout())

### Name: writeLD
### Title: Linkage Disequilibrium (LD) computation with PLINK
### Aliases: writeLD
### Keywords: Enrichment initFiles LD ld linkage Disequilibrium PLINK

### ** Examples

## Not run:
# signalFile <- "./extData/Signal/toySignal.txt"
# snpInfoDir <- "./extData/snpInfo/"
# snpListDir <- "./extData/List/"
# writeLD(pattern = "Chrom", snpInfoDir, snpListDir,
#         signalFile, ldThresh = 0.8, onlySignal = TRUE, 
#         mc.cores = detectCores())
## End (Not run)



### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
