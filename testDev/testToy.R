rm(list = ls())
# library(snpEnrichment)
source("/ep10/disks/SANA8/Mickael/ENRICHMENT/Scripts/main.R")
setwd("/ep10/disks/SANA8/Mickael/ENRICHMENT/Run/")


### Example
# path <- path.package("snpEnrichment")
# snpListDir <- paste(path, "/extdata/List/", sep = "")
# signalFile <- paste(path, "/extdata/Signal/toySignal.txt", sep = "")
# excludeFile <- paste(path, "/extdata/Exclude/toyExclude.txt", sep = "")
# snpInfoDir <- paste(path, "/extdata/snpInfo/", sep = "")
# initFiles(pattern = "Chrom", snpInfoDir, signalFile, ldThresh = 0.08, LD = TRUE)

pattern <- "Chrom"

transcriptFile <- "/ep10/disks/SANA8/Mickael/ENRICHMENT/Data/Toy/transcript.txt"
# transcriptFile <- FALSE
# data(transcript)
# transcriptFile <- transcript

sigThresh <- 0.05
ldThresh <- 0.8
nSample <- 1000
MAFpool <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5)
distThresh <- 1000
LD <- TRUE
extendMethod <- "ld"

snpListDir <- "/ep10/disks/SANA8/Mickael/ENRICHMENT/Data/Toy/List/"
signalFile <- "/ep10/disks/SANA8/Mickael/ENRICHMENT/Data/Toy/Signal/toySignal.txt"
excludeFile <- "/ep10/disks/SANA8/Mickael/ENRICHMENT/Data/Toy/Exclude/toyExclude.txt"
snpInfoDir <- "/ep10/disks/SANA8/Mickael/ENRICHMENT/Data/Toy/snpInfo/"

mc.cores <- detectCores()

initFiles(pattern, snpInfoDir, signalFile, ldThresh, LD, mc.cores)
toy_M1 <- readEnrichment(pattern, signalFile, transcriptFile, snpListDir, snpInfoDir, distThresh, sigThresh, LD, extendMethod, mc.cores)
toyR_M1 <- toy_M1
# save(toyR_M1, file = "toyR_M1.RData")
# toyM1 <- toy_M1
# save(toyM1, file = "toyM1.RData")
# load("testRead_M1.RData")

reSample(object = toy_M1, nSample, sigThresh, MAFpool, extendMethod, mc.cores)
toyC_M1 <- toy_M1
save(toyC_M1, file = "toyC_M1.RData")
# load("testCompute_M1.RData")

toyE_M1 <- excludeSNP(toy_M1, excludeFile, nSample, sigThresh, MAFpool, extendMethod, mc.cores)
save(toyE_M1, file = "toyE_M1.RData")
# load("testExclude_M1.RData")

result1 <- compareEnrichment(toyC_M1, toyE_M1, pattern, signalFile, transcriptFile, snpInfoDir, distThresh, LD, ldThresh, nSample, sigThresh, MAFpool, extendMethod, mc.cores)


extendMethod = "block"
toy_M2 <- readEnrichment(pattern, signalFile, transcriptFile, snpListDir, snpInfoDir, distThresh, sigThresh, LD, extendMethod, mc.cores)
toyR_M2 <- toy_M2
# save(toyR_M2, file = "toyR_M2.RData")
# toyM2 <- toy_M2
# save(toyM2, file = "toyM2.RData")
# load("testRead_M1.RData")

reSample(object = toy_M2, nSample, sigThresh, MAFpool, extendMethod, mc.cores)
toyC_M2 <- toy_M2
save(toyC_M2, file = "toyC_M2.RData")
# load("testCompute_M1.RData")

toyE_M2 <- excludeSNP(toy_M2, excludeFile, nSample, sigThresh, MAFpool, extendMethod, mc.cores)
save(toyE_M2, file = "toyE_M2.RData")
# load("testExclude_M1.RData")

result2 <- compareEnrichment(toyC_M2, toyE_M2, pattern, signalFile, transcriptFile, snpInfoDir, distThresh, LD, ldThresh, nSample, sigThresh, MAFpool, extendMethod, mc.cores)

