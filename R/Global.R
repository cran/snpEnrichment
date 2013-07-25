setGeneric(name = "newSNP", def = function(List, Table, EnrichmentRatio, Z, PValue, Resampling){standardGeneric("newSNP")})
setGeneric(name = "chromosome", def = function(Data, LD, eSNP, xSNP){standardGeneric("chromosome")})
setGeneric(name = "enrichment", def = function(Loss, Parameters, eSNP, xSNP, Chromosomes){standardGeneric("enrichment")})
setGeneric(name = "doLDblock", def = function(object, mc.cores = detectCores()){standardGeneric("doLDblock")})
setGeneric(name = "excludeSNP", def = function(object, excludeFile, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()){standardGeneric("excludeSNP")})
setGeneric(name = "computeER", def = function(object, sigThresh = 0.05, mc.cores = detectCores()){standardGeneric("computeER")})
setGeneric(name = "reSample", def = function(object, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()){standardGeneric("reSample")})
setGeneric(name = "compareEnrichment", def = function(list1, list2, pattern = "Chrom", signalFile, transcriptFile = FALSE, snpInfoDir, distThresh = 1000, LD = FALSE, ldThresh = 0.8, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()){standardGeneric("compareEnrichment")})
setGeneric(name = "is.enrichment", def = function(object){standardGeneric("is.enrichment")})
setGeneric(name = "is.chromosome", def = function(object){standardGeneric("is.chromosome")})
setGeneric(name = "is.SNP", def = function(object){standardGeneric("is.SNP")})
setGeneric(name = "reset", def = function(object, i){standardGeneric("reset")})
# setGeneric(name = "summary", def = function(object, extended = TRUE, complete = TRUE){standardGeneric("summary")})

setMethod(f = "reSample", signature = "ANY", definition = function(object, nSample, sigThresh, MAFpool, extendMethod, mc.cores) {
    if (!(is.enrichment(object) & is.chromosome(object))){
        stop('[Method:reSample] not available for "', class(object), '" object', call. = FALSE)
        return(invisible())
    } else {}
})
setMethod(f = "excludeSNP", signature = "ANY", definition = function(object, excludeFile, extendMethod, mc.cores) {
    if (!is.enrichment(object)){
        stop('[Method:excludeSNP] not available for "', class(object), '" object', call. = FALSE)
        return(invisible())
    } else {}
})
setMethod(f = "reset", signature = "ANY", definition = function(object, i) {
    if (!(is.enrichment(object) & is.chromosome(object))){
        stop('[Method:reset] not available for "', class(object), '" object', call. = FALSE)
        return(invisible())
    } else {}
})


.verbose <- function(expr) {invisible(capture.output(expr))}


GC <- function(verbose = getOption("verbose"), reset = FALSE) {
    while (!identical(gc(verbose, reset)[, 4], gc(verbose, reset)[, 4])) {}
    return(gc(verbose, reset))
}


.checkFilePath <- function(path) {
    END <- unlist(regmatches(path, regexec(".$", path)))
    # START <- unlist(regmatches(path, regexec("^.", path)))
    # if (START != "/") {
        # path <- paste("/", path, sep = "")
    # } else {}
    if (END != "/") {
        path <- paste(path, "/", sep = "")
    } else {}
    return(path)
}


mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, 
                      mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
    if (Sys.info()[["sysname"]] != "Windows") {
        mc.cores.old <- mc.cores
        sysMemFree <- system("egrep 'MemFree' /proc/meminfo", intern = TRUE)
        sysMemAvailable <- 0.90*as.numeric(unlist(regmatches(sysMemFree, regexec("[0-9]+", sysMemFree))))
        sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"))[8])
        mc.cores <- min(as.integer(mc.cores), floor(sysMemAvailable/sysProc))
        if (mc.cores != mc.cores.old) {
            msg <- paste('[mclapply2] To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".', sep = "")
            warning(msg, call. = FALSE)
        } else {}
    } else {}
    require(parallel)
    return(mclapply(X = X, FUN = FUN, ..., 
                    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent, 
                    mc.cores = mc.cores, mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive))
}


maxCores <- function(mc.cores = detectCores()) {
    sysMemTmp <- system("head -1 /proc/meminfo", intern = TRUE)
    sysMem <- as.numeric(unlist(regmatches(sysMemTmp, regexec("[0-9]+", sysMemTmp))))
    sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"))[8])
    return(min(mc.cores, floor(sysMem/sysProc)))
}


.showArgs <- function(args) {
    res <- NULL
    for (iArg in seq(length(args[[1]]))) {
        res <- c(res, paste(names(args[[1]][iArg]), args[[1]][iArg], sep=" = "))
    }
    paste(names(args), '(', paste(res, collapse=", "), ')', sep="")
}


# .writeSignal <- function(pattern, snpInfoDir, snpListDir, signalFile) {
.writeSignal <- function(pattern, snpInfoDir, signalFile) {
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t")))>1) {
        signal <- read.delim(file = signalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                        colClasses = c("character", "numeric"), na.string = c("NA", ""), 
                        check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
    } else {
        if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " ")))>1) {
            signal <- read.delim(file = signalFile, header = TRUE, sep = " ", stringsAsFactors = FALSE, 
                        colClasses = c("character", "numeric"), na.string = c("NA", ""), 
                        check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
        } else {
            stop('[Enrichment:initFiles] only " " and "\t" are allowed as columns separator in Signal file', call. = FALSE)
            return(invisible())
        }
    }

    chrom.bim <- read.delim(file = paste(snpInfoDir, pattern, ".bim", sep = ""), header = FALSE, stringsAsFactors = FALSE, 
                            colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""), 
                            check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "", "", "", ""))
    signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
    # eSNP <- .readSNP(pattern = pattern, snpListDir = snpListDir)
    # eSNPplink <- matrix(unique(merge(eSNP, chrom.bim, by = "SNP")[, 1]), dimnames = list(NULL, "SNP"))
    # signal <- merge(signal, eSNPplink, by = "SNP", all = TRUE)[, c(2, 1, 3)]
    eval(parse(text = paste('write.table(signal, file = "', 
                            paste(snpInfoDir, pattern, ".signal", sep = ""), 
                            '", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")', sep = "")))
    return(invisible())
}


writeLD <- function(pattern = "Chrom", snpInfoDir, signalFile, ldThresh = 0.8, onlySignal = TRUE, mc.cores = detectCores()) {
    if (missing(pattern) | missing(snpInfoDir) | missing(signalFile) | missing(ldThresh)) {
        stop('[Enrichment:writeLD] argument(s) missing', call. = FALSE)
        return(invisible())
    } else {}
    snpInfoDir <- .checkFilePath(snpInfoDir)
    if (Sys.info()[["sysname"]] == "Linux") {
        if (as.logical(system("rpm -qa | grep plink > /dev/null && echo TRUE || echo FALSE", intern = TRUE))) {
            if (onlySignal) {
                FILES <- list.files(snpInfoDir)
                resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
                    newPattern <- unlist(strsplit(grep(paste(pattern, iChr, ".*.bim", sep = ""), FILES, value = TRUE), ".bim"))[1]
                    .writeSignal(pattern = newPattern, snpInfoDir, signalFile)
                    system(paste("dir=", eval(paste(snpInfoDir, newPattern, sep = "")), 
                                "\nsnplist=", eval(paste(snpInfoDir, newPattern, ".signal", sep = "")), 
                                "\nldThresh=", eval(ldThresh), 
                                "\nplink --silent --noweb --bfile $dir --extract $snplist --r2 --ld-snp-list $snplist --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 $ldThresh --out $dir \n", sep = ""))
                })
            } else {
                FILES <- list.files(snpInfoDir)
                resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
                    newPattern <- unlist(strsplit(grep(paste(pattern, iChr, ".*.bim", sep = ""), FILES, value = TRUE), ".bim"))[1]
                    .writeSignal(pattern = newPattern, snpInfoDir, signalFile)
                    system(paste("dir=", eval(paste(snpInfoDir, newPattern, sep = "")), 
                                "\nsnplist=", eval(paste(snpInfoDir, newPattern, ".signal", sep = "")), 
                                "\nldThresh=", eval(ldThresh), 
                                "\nplink --silent --noweb --bfile $dir --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 $ldThresh --out $dir \n", sep = ""))
                })
            }
            return(invisible())
        } else {
            stop('[Enrichment:initFiles] "PLINK" is not installed', call. = FALSE)
            return(invisible())
        }
    } else {
        stop('[Enrichment:initFiles] only available on "UNIX" system', call. = FALSE)
        return(invisible())
    }
}


.writeFreq <- function(pattern, snpInfoDir) {
    if (Sys.info()[["sysname"]] == "Linux") {
        if (as.logical(system("rpm -qa | grep plink > /dev/null && echo TRUE || echo FALSE", intern = TRUE))) {
            system(paste("dir=", eval(paste(snpInfoDir, pattern, sep = "")), "\nplink --silent --noweb --bfile $dir --freq --out $dir \n", sep = ""))
            freq <- read.delim(file = paste(snpInfoDir, pattern, ".frq", sep = ""), header = TRUE, stringsAsFactors = FALSE, 
                                colClasses = c("NULL", "character", "NULL", "NULL", "numeric", "NULL"), na.string = c("NA", ""), 
                                check.names = FALSE, strip.white = TRUE, sep = "")
            bim <- read.delim(file = paste(snpInfoDir, pattern, ".bim", sep = ""), header = FALSE, stringsAsFactors = FALSE, 
                                colClasses = c("numeric", "character", "NULL", "numeric", "NULL", "NULL"), na.string = c("NA", ""), 
                                check.names = FALSE, strip.white = TRUE, col.names = c("CHR", "SNP", "", "POS", "", ""))
            plinkData <- merge(freq, bim, by = "SNP")[, c("CHR", "SNP", "POS", "MAF")]
            write.table(plinkData, paste(snpInfoDir, pattern, ".all", sep = ""), row.names = FALSE, sep = "\t")
            return(invisible())
        } else {
            stop('[Enrichment:initFiles] "PLINK" is not installed', call. = FALSE)
            return(invisible())
        }
    } else {
        stop('[Enrichment:initFiles] only available on "UNIX" system. \nPlease run "plink --freq" for all chromosomes.', call. = FALSE)
        return(invisible())
    }
}


initFiles <- function(pattern = "Chrom", snpInfoDir, snpListDir, signalFile, ldThresh = 0.8, LD = FALSE, mc.cores = detectCores()) {
    if (missing(snpInfoDir) | missing(signalFile) | missing(ldThresh)) {
        stop('[Enrichment:initFiles] argument(s) missing', call. = FALSE)
        return(invisible())
    } else {}
    snpInfoDir <- .checkFilePath(snpInfoDir)
    FILES <- list.files(snpInfoDir)
    cat("All files are ready for chromosome:\n  ")
    resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(jChr){
        newPattern <- unlist(strsplit(grep(paste(pattern, jChr, ".*.bim", sep = ""), FILES, value = TRUE), ".bim"))[1]
        .writeSignal(pattern = newPattern, snpInfoDir = snpInfoDir, signalFile = signalFile)
        .writeFreq(pattern = newPattern, snpInfoDir = snpInfoDir)
        if (LD) {
            if (Sys.info()[["sysname"]] == "Linux") {
                if (as.logical(system("rpm -qa | grep plink > /dev/null && echo TRUE || echo FALSE", intern = TRUE))) {
                    res2 <- system(paste("dir=", eval(paste(snpInfoDir, newPattern, sep = "")), 
                                        "\nsnplist=", eval(paste(snpInfoDir, newPattern, ".signal", sep = "")), 
                                        "\nldThresh=", eval(ldThresh), 
                                        "\nplink --silent --noweb --bfile $dir --extract $snplist --r2 --ld-snp-list $snplist --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 $ldThresh --out $dir \n", sep = ""))
                } else {
                    stop('[Enrichment:initFiles] "PLINK" is not installed', call. = FALSE)
                    return(invisible())
                }
            } else {
                stop('[Enrichment:initFiles] only available on "UNIX" system', call. = FALSE)
                return(invisible())
            }
        } else {}
        cat(jChr, " ", sep = "")
        return(invisible())
    })
    cat("\n")
    return(invisible())
}


.readSignal <- function(pattern, snpInfoDir) {
    signal <- read.delim(file = paste(snpInfoDir, pattern, ".signal", sep = ""), header = TRUE, stringsAsFactors = FALSE, 
                            colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""), 
                            check.names = FALSE, strip.white = TRUE, col.names =  c("", "SNP", "PVALUE"))
    return(signal)
}


.readSNP <- function(pattern, snpListDir) {
    snpListFile <- paste(snpListDir, grep(pattern, list.files(snpListDir), value = TRUE), sep = "")[1]
    snpList <- read.delim(file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
                            na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
    switch(EXPR = as.character(ncol(snpList)),
        "1" = {colnames(snpList) <- "SNP"},
        "2" = {colnames(snpList) <- c("CHR", "SNP")},
        "3" = {colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT")},
        colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
    )
    if (ncol(snpList)>=3) {
        snpList <- snpList[, c("SNP", "TRANSCRIPT")]
    } else {}
    if (nrow(snpList)!=0) {
        snpList[, "eSNP"] <- 1
    } else {}
    return(snpList)
}


.readTranscript <- function(transcriptFile) {
    if (all(class(try(close(file(transcriptFile)), silent = TRUE))!="try-error")) {
        transcript <- read.delim(file = transcriptFile, header = TRUE, stringsAsFactors = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
        transcript <- transcript[, 1:4]
        colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
    } else {
        if (class(transcriptFile) %in% c("matrix", "data.frame")) {
            transcript <- transcriptFile[, 1:4]
            colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
        } else {
            stop('[Enrichment:readEnrichment] "transcriptFile" required "matrix", "data.frame" or "txt file"', call. = FALSE)
            return(invisible())
        }
    }
    return(transcript)
}


.readFreq <- function(pattern, snpInfoDir) {
    freq <- read.delim(file = paste(snpInfoDir, pattern, ".all", sep = ""), header = TRUE, stringsAsFactors = FALSE, 
                        colClasses = c("numeric", "character", "numeric", "numeric"), na.string = c("NA", ""), 
                        check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "POS", "MAF"))
    return(freq)
}


.readLD <- function(pattern, snpInfoDir) {
    ldData <- read.delim(file = paste(snpInfoDir, pattern, ".ld", sep = ""), header = TRUE, stringsAsFactors = FALSE, 
                            colClasses = c("NULL", "NULL", "character", "NULL", "numeric", "character", "NULL"), na.string = c("NA", ""), 
                            check.names = FALSE, strip.white = TRUE, sep = "")
    return(ldData)
}


.readFiles <- function(pattern, snpInfoDir, snpListDir, transcriptFile, distThresh, LD) {
    FILES <- list.files(snpInfoDir)
    newPattern <- unlist(strsplit(grep(paste(pattern, ".*.bim", sep = ""), FILES, value = TRUE), ".bim"))[1]
    eSNP <- .readSNP(pattern = newPattern, snpListDir = snpListDir)
    signal <- .readSignal(pattern = newPattern, snpInfoDir = snpInfoDir)
    plinkData <- .readFreq(pattern = newPattern, snpInfoDir = snpInfoDir)
    
    if (LD) {
        linkageData <- .readLD(pattern = newPattern, snpInfoDir = snpInfoDir)
    } else {
        linkageData <- character()
    }

    signalPlink <- merge(signal, plinkData, by = "SNP")
    if (nrow(eSNP)!=0) {
        eSNPunique <- eSNP[!duplicated(eSNP[, "SNP"]), ]
        snpSignal <- merge(signalPlink, eSNPunique[, c("SNP", "eSNP")], by = "SNP", all.x = TRUE)
        snpSignal[, "eSNP"][is.na(snpSignal[, "eSNP"])] <- 0
    } else {
        snpSignal <- signalPlink
        snpSignal[, "eSNP"] <- 0
    }
    
    chrNumber <- unlist(regmatches(newPattern, regexec("[0-9]+", newPattern)))
    if (any(transcriptFile != FALSE)) {
        transcript <- .readTranscript(transcriptFile = transcriptFile)
        transcriptCHR <- na.exclude(transcript[transcript[, "CHR"] == unique(snpSignal[, "CHR"]), c("START", "END")])

        cisFunc <- function(line, distThresh) {
            position <- as.numeric(line[4])
            CIS <- any(position > (transcriptCHR[, "START"] - distThresh) & (position < (transcriptCHR[, "END"] + distThresh)))
            return(c(line, CIS = as.numeric(CIS)))
        }

        temp <- as.data.frame(t(apply(snpSignal, 1, cisFunc, distThresh = distThresh*10^3)), stringsAsFactors = FALSE)
        temp <- unique(temp[, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP", "CIS")])

        data <- transform(temp, SNP = as.character(temp$SNP), PVALUE = as.numeric(temp$PVALUE), CHR = as.numeric(temp$CHR), POS = as.numeric(temp$POS), MAF = as.numeric(temp$MAF), eSNP = as.numeric(temp$eSNP), CIS = as.numeric(temp$CIS))
        signalLoss <- c(NA, NA, length(unique(data[, "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])), sum(data[, "CIS"]))
        data <- data[data[, "CIS"] == 1, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")]
    } else {
        temp <- unique(snpSignal[, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")])
        data <- transform(temp, SNP = as.character(temp$SNP), PVALUE = as.numeric(temp$PVALUE), CHR = as.numeric(temp$CHR), POS = as.numeric(temp$POS), MAF = as.numeric(temp$MAF), eSNP = as.numeric(temp$eSNP))
        signalLoss <- c(NA, NA, length(unique(data[, "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])))
    }
    data[, "xSNP"] <- 0
    if (any(duplicated(data[, "SNP"]))) {
        cat(data[, "SNP"][duplicated(data[, "SNP"])], "\n")
        stop('[Enrichment:readEnrichment] Duplicated SNPs in Signal', call. = FALSE)
        return(invisible())
    } else {
        rownames(data) <- data[, "SNP"]
    }
    snpLoss <- c(length(eSNP[, "SNP"]), length(unique(eSNP[, "SNP"])), length(unique(intersect(eSNP[, "SNP"], signal[, "SNP"]))), length(unique(intersect(eSNP[, "SNP"], signal[!is.na(signal[, "PVALUE"]), "SNP"]))), length(data[data[, "eSNP"]==1 & !is.na(data[, "PVALUE"]), "SNP"]))
    return(list(data = data, LD = linkageData, snpLoss = snpLoss, signalLoss = signalLoss))
}


readEnrichment <- function(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 1000, sigThresh = 0.05, LD = FALSE, extendMethod = "ld", mc.cores = detectCores()) {
    cat("############# Read Enrichment ##############\n")
    if (missing(signalFile) | missing(snpListDir) | missing(snpInfoDir)) {
        stop('[Enrichment:readEnrichment] argument(s) missing', call. = FALSE)
        return(invisible())
    } else {}
    snpInfoDir <- .checkFilePath(snpInfoDir)
    snpListDir <- .checkFilePath(snpListDir)
    cat("  Read Chromosomes:\n    ")
    if (extendMethod == "block") {
        cat("- do LD block ...\n        ")
    } else {}
    resParallel <- mclapply2(seq(22), mc.cores = min(22, mc.cores), FUN = function(jChr){
        files <- .readFiles(pattern = paste(pattern, jChr, sep = ""), snpInfoDir = snpInfoDir, snpListDir = snpListDir, transcriptFile = transcriptFile, distThresh = distThresh, LD = LD)
        if (LD) {
            linkageData <- files$LD
            dataTmp <- files$data
            data <- dataTmp[order(dataTmp$POS), ]
            snpSignal <- data[, "SNP"]
            
            linkageData <- linkageData[linkageData[, "SNP_A"]%in%snpSignal & linkageData[, "SNP_B"]%in%snpSignal, ]
            linkageData <- linkageData[order(linkageData[, "BP_B"]), ]
            ldData <- c(linkageData[, 1], linkageData[, 3])
            names(ldData) <- c(linkageData[, 3], linkageData[, 1])
            if (extendMethod == "ld") {
                eSNP <- data[data[, "eSNP"]==1, "SNP"]
                xSNP <- intersect(data[, "SNP"], unique(c(data[data[, "eSNP"]==1, "SNP"], ldData[names(ldData) %in% eSNP])))
                data[intersect(xSNP, snpSignal), "xSNP"] <- 1
                data <- data[!is.na(data$PVALUE), ]
                resChr <- chromosome(Data = data, LD = ldData)
            } else {
                filesTmp <- files[3:4]
                rm(files, dataTmp, linkageData, snpSignal)
                GC()
                if (extendMethod == "block") { # extended LD for SNP which are not signal!!!!
                    resTmp <- chromosome(Data = data, LD = ldData)
                    resTmp <- doLDblock(object = resTmp, mc.cores = mc.cores)
                    data <- resTmp@Data
                    # rm(resTmp)
                    # GC()
                    data[data[, "IDBLOCK"]%in%unique(data[data[, "eSNP"]==1, "IDBLOCK"]), "xSNP"] <- 1
                    data <- data[!is.na(data[, "PVALUE"]), ]
                    files <- filesTmp
                    resChr <- chromosome(Data = data, LD = ldData)
                } else {
                    stop('[Enrichment:readEnrichment] wrong method for extend eSNP list', call. = FALSE)
                    return(invisible())
                }
            }
        } else {
            resChr <- chromosome(Data = files$data, LD = files$LD)
        }
        cat(jChr, " ", sep = "")
        return(list(resChr, files$snpLoss, files$signalLoss))
    })
    names(resParallel) <- paste("Chrom", seq(22), sep="")
    cat("\n")
    
    result = enrichment()
    result@Chromosomes <- lapply(resParallel, function(l) {return(l[[1]])})
    snpLoss <- t(as.data.frame(lapply(resParallel, function(l) {return(l[[2]])})))
    signalLoss <- apply(t(as.data.frame(lapply(resParallel, function(l) {return(l[[3]])}))), 2, sum)
    loss <- rbind(Signal = signalLoss, AllChrom = apply(snpLoss, 2, sum), snpLoss)
    colnames(loss) <- c("Rows", "Unique", "Intersect.Ref", "Intersect.Signal", "CIS")
    rm(resParallel)
    
    SNPs <- result["List", seq(22)]
    result@eSNP@List <- SNPs[["eSNP"]]
    result@xSNP@List <- SNPs[["xSNP"]]
    rm(SNPs)

    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t")))>1) {
        signal <- read.delim(file = signalFile, header = TRUE, sep = "\t", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"), stringsAsFactors = FALSE)
    } else {
        if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " ")))>1) {
            signal <- read.delim(file = signalFile, header = TRUE, sep = " ", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"), stringsAsFactors = FALSE)
        } else {
            stop('[Enrichment:readEnrichment] only " " and "\t" are allowed as columns separator in Signal file', call. = FALSE)
            return(invisible())
        }
    }
    loss[1, c(1, 2)] <- c(length(signal[, "SNP"]), length(unique(signal[, "SNP"])))
    result@Loss <- as.data.frame(loss)
    rm(signal)

    sysCall <- sys.call()
    argsSNP <- as.list(sysCall[-1])
    formal <- as.list(names(formals(as.character(sysCall))))
    names(formal) <- formal
    for (iArg in names(formal)) {
        if (identical(grep(iArg, names(argsSNP)), integer(0))) {
            if (grep(iArg, names(formal)) > length(argsSNP)) {
                formal[[iArg]] <- eval(parse(text = formal[[iArg]]))
            } else {
                argTmp <- argsSNP[[grep(iArg, names(formal))]]
                if (is.character(argTmp) | is.matrix(argTmp) | is.data.frame(argTmp)) { # check if transcript is still evaluated
                    formal[[iArg]] <- argsSNP[[grep(iArg, names(formal))]]
                } else {
                    formal[[iArg]] <- eval(parse(text = argsSNP[[grep(iArg, names(formal))]]))
                }
            }
        } else {
            formal[[iArg]] <- argsSNP[[iArg]]
        }
    }
    result@Parameters$read <- formal
    
    result <- computeER(object = result, sigThresh = sigThresh, mc.cores = mc.cores)
    cat("########### Read Enrichment Done ###########\n")
    return(result)
}


.enrichmentRatio <- function(table) {
    tmp <- as.numeric(table)
    return((tmp[1]*tmp[4])/(tmp[2]*tmp[3]))
}

.reSample <- function(object, nSample, sigThresh, MAFpool, extendMethod, mc.cores) {
    nResampling <- nrow(object@eSNP@Resampling)
    data <- object["Data"]
    chrLD <- object["LD"]
    isLD <- length(chrLD) != 0
    eEnrichment <- object@eSNP@Table
    eEnrichRatio <- object@eSNP@EnrichmentRatio
    eSNPlist <- object@eSNP@List
    data[, "MAFpool"] <- NA
    data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
    nPool <- nlevels(data[, "MAFpool"])
    if (extendMethod == "block") {
        data[, "BLOCKpool"] <- NA
        data[, "BLOCKpool"] <- as.factor(cut(data[, "LENGTH"], breaks = unique(c(0, 1, floor(mean(sort(unique(data[, "LENGTH"]))[-1])), max(data[, "LENGTH"]))), labels = FALSE, include.lowest = TRUE))
        eSNPlistPool <- table(data[data[, "eSNP"] == 1, "BLOCKpool"], data[data[, "eSNP"] == 1, "MAFpool"])
    } else {
         eSNPlistPool <- table(data[data[, "eSNP"] == 1, "MAFpool"])
    }
    
    if (nResampling == 0) {
        eSNPsample <- NULL
        nSampleMin <- min(nSample, max(1000, nSample*10/100))
    } else {
        eSNPsample <- object@eSNP@Resampling
        nSampleMin <- nSample
    }
    
    if (isLD) {
        if (nResampling == 0) {xSNPsample <- NULL} else {xSNPsample <- object@xSNP@Resampling}
        xEnrichment <- object@xSNP@Table
        xEnrichRatio <- object@xSNP@EnrichmentRatio
        xSNPlist <- object@xSNP@List
        if (extendMethod == "block") {
            xSNPlistPool <- table(data[data[, "xSNP"] == 1, "BLOCKpool"], data[data[, "xSNP"] == 1, "MAFpool"], deparse.level = 0)
            popSNP4Sample <- lapply(split(data[data[, "xSNP"] != 1, c("SNP", "BLOCKpool")], data[data[, "xSNP"] != 1, "MAFpool"]), function(mafList){
                split(mafList[, "SNP"], mafList[, "BLOCKpool"])
            })
        } else {
            xSNPlistPool <- table(data[data[, "xSNP"] == 1, "MAFpool"], deparse.level = 0)
            popSNP4Sample <- split(data[data[, "xSNP"] != 1, "SNP"], data[data[, "xSNP"] != 1, "MAFpool"])
        }
    } else {
        popSNP4Sample <- split(data[data[, "eSNP"] != 1, "SNP"], data[data[, "eSNP"] != 1, "MAFpool"])
    }
    assoc <- factor(data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    SNPlist <- data[, "SNP"]
    if (extendMethod == "block") {
        data0 <- data[, c("SNP", "IDBLOCK")]
    } else {}
    rm(data)
    GC()
    
    cat(0, ".. ", sep = "")
    resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function(i){
        if (extendMethod == "block") {
            eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){
                lapply(seq(length(popSNP4Sample[[g]])), function(b){
                    sample(popSNP4Sample[[g]][[b]], min(eSNPlistPool[b, g], length(popSNP4Sample[[g]][[b]])))
                })
            }))
        } else {
            eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
        }
        eSNPenrichStats <- table(assoc, factor(SNPlist%in%eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
        if (isLD) {
            if (extendMethod == "ld") {
                xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
            } else {
                xSNPlistRandom <- data0[data0[, "IDBLOCK"] %in% data0[data0[, "SNP"] %in% eSNPlistRandom, "IDBLOCK"], "SNP"]
            }
            xSNPenrichStats <- table(assoc, factor(SNPlist%in%xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
            xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
        } else {
            xTmp <- NULL
        }
        eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
        return(list(eSNP = eTmp, xSNP = xTmp))
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l){l$eSNP})))
    eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
    Ze <- (eEnrichRatio-eMeanEnrichRatio)/((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/length(eSNPsample[, 5]))
    # Ze <- (eEnrichRatio-mean(eSNPsample[, 5]))/((sum((eSNPsample[, 5]-mean(eSNPsample[, 5]))^2))/length(eSNPsample[, 5]))
    if (isLD) {
        xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP})))
        xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
        Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
        # Zx <- (xEnrichRatio-mean(xSNPsample[, 5]))/((sum((xSNPsample[, 5]-mean(xSNPsample[, 5]))^2))/length(xSNPsample[, 5]))
    } else {
        Zx <- 1
    }
    iSample <- nSampleMin
    cat(iSample, ".. ", sep = "")
    catStep <- (nSample/10)
    rm(resParallel)
    GC()
    
    while (iSample<nSample & ((Ze>0 & Ze<15) | is.na(Ze) | (Zx>0 & Zx<15) | is.na(Zx))) {
        resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function(i){
            if (extendMethod == "block") {
                eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){
                    lapply(seq(length(popSNP4Sample[[g]])), function(b){
                        sample(popSNP4Sample[[g]][[b]], min(eSNPlistPool[b, g], length(popSNP4Sample[[g]][[b]])))
                    })
                }))
            } else {
                eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
            }
            eSNPenrichStats <- table(assoc, factor(SNPlist%in%eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
            if (isLD) {
                if (extendMethod == "ld") {
                    xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
                } else {
                    xSNPlistRandom <- data0[data0[, "IDBLOCK"] %in% data0[data0[, "SNP"] %in% eSNPlistRandom, "IDBLOCK"], "SNP"]
                }
                xSNPenrichStats <- table(assoc, factor(SNPlist%in%xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
                xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
            } else {
                xTmp <- NULL
            }
            eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
            return(list(eSNP = eTmp, xSNP = xTmp))
        })
        eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l){l$eSNP})))
        eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
        Ze <- (eEnrichRatio-eMeanEnrichRatio)/((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/length(eSNPsample[, 5]))
        # Ze <- (eEnrichRatio-mean(eSNPsample[, 5]))/((sum((eSNPsample[, 5]-mean(eSNPsample[, 5]))^2))/length(eSNPsample[, 5]))
        if (isLD) {
            xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP})))
            xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
            Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
            # Zx <- (xEnrichRatio-mean(xSNPsample[, 5]))/((sum((xSNPsample[, 5]-mean(xSNPsample[, 5]))^2))/length(xSNPsample[, 5]))
        } else {
            Zx <- 1
        }
        iSample <- iSample + nSampleMin
        if ((iSample-(iSample%/%catStep*catStep)) == 0) {cat(iSample, ".. ", sep = "")}
        rm(resParallel)
    }
    # cat("\n")


    eSNP <- newSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = pnorm(Ze, lower.tail = FALSE), Resampling = eSNPsample)
    if (isLD) {
        xSNP <- newSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = pnorm(Zx, lower.tail = FALSE), Resampling = xSNPsample)
    } else {
        xSNP <- newSNP()
    }
    object@eSNP <- eSNP
    object@xSNP <- xSNP
    return(object)
}


.compareEnrich <- function(object1, object2, nSample, sigThresh, MAFpool, extendMethod, mc.cores) {
    DATA <- object1["Data"]
    chrLD <- object1["LD"]
    isLD <- length(chrLD) != 0
    eEnrichment <- object1@eSNP@Table
    eEnrichRatio <- object1@eSNP@EnrichmentRatio
    eSNPlist <- object1@eSNP@List
    eList <- union(object1@eSNP@List, object2@eSNP@List) # All eSNP in object1 and object2
    if (isLD) {
        xList <- union(object1@xSNP@List, object2@xSNP@List)
        data <- DATA[DATA[, "SNP"] %in% union(eList, xList), ]
    } else {
        data <- DATA[DATA[, "SNP"] %in% eList, ]
    }
    
    data[, "MAFpool"] <- NA
    data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
    nPool <- nlevels(data$MAFpool)
    eSNPlistPool <- table(data[eSNPlist, "MAFpool"])

    eSNPsample <- NULL
    nSampleMin <- nSample
    if (isLD) {
        xEnrichment <- object1@xSNP@Table
        xEnrichRatio <- object1@xSNP@EnrichmentRatio
        xSNPlist <- object1@xSNP@List
        xSNPlistPool <- table(data[xSNPlist, "MAFpool"])
        xSNPsample <- NULL
    } else {}
    
    # SNPlist <- data[, "SNP"]
    popSNP4Sample <- split(eList, data[eList, "MAFpool"])
    assoc_eSNP <- factor(data[eList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    assoc_xSNP <- factor(data[xList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    if (extendMethod == "block") {
        data0 <- data[, c("SNP", "IDBLOCK")]
    } else {}
    rm(data)
    
    cat("0.. ")
    resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function(i){
        eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
        eSNPenrichStats <- table(assoc_eSNP, factor(eList%in%eSNPlistRandom, levels = c(FALSE, TRUE)))
        eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
        if (isLD) {
            if (extendMethod == "ld") {
                xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
            } else {
                xSNPlistRandom <- data0[data0[, "IDBLOCK"] %in% data0[data0[, "SNP"] %in% eSNPlistRandom, "IDBLOCK"], "SNP"]
            }
            xSNPenrichStats <- table(assoc_xSNP, factor(xList%in%xSNPlistRandom, levels = c(FALSE, TRUE)))
            xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
        } else {
            xTmp <- NULL
        }
        return(list(eSNP = eTmp, xSNP = xTmp))
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l){l$eSNP})))

    eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
    if (is.matrix(eSNPanyDup)) {
        if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
            eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
        } else {}
        if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
            eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[!is.infinite(eSNPsample[, 5]), 5])*1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
        } else {}
        eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
        Ze <- (eEnrichRatio-eMeanEnrichRatio)/((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/length(eSNPsample[, 5]))
        # Ze <- (eEnrichRatio-mean(eSNPsample[, 5]))/((sum((eSNPsample[, 5]-mean(eSNPsample[, 5]))^2))/length(eSNPsample[, 5]))
    } else {
        Ze <- as.numeric(NA)
    }
    
    if (isLD) {
        xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP})))
        xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
        if (is.matrix(xSNPanyDup)) {
            if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
                xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
            } else {}
            if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
                xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[!is.infinite(xSNPsample[, 5]), 5])*1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
            } else {}
            xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
            Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
            # Zx <- (xEnrichRatio-mean(xSNPsample[, 5]))/((sum((xSNPsample[, 5]-mean(xSNPsample[, 5]))^2))/length(xSNPsample[, 5]))
        } else {
            Zx <- as.numeric(NA)
        }
    } else {
        Zx <- as.numeric(NA)
    }
    iSample <- nSampleMin
    cat(iSample, ".. ", sep = "")
    catStep <- (nSample/10)
    rm(resParallel)
    GC()
    
    while (iSample<nSample & ((Ze>0 & Ze<15) | is.na(Ze) | (Zx>0 & Zx<15) | is.na(Zx))) {
        resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function(i){
            eSNPlistRandom <- unlist(sapply(seq(nPool), function(g){sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
            eSNPenrichStats <- table(assoc_eSNP, factor(eSNPlist%in%eSNPlistRandom, levels = c(FALSE, TRUE)))
            eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
            if (isLD) {
                if (extendMethod == "ld") {
                    xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
                } else {
                    xSNPlistRandom <- data0[data0[, "IDBLOCK"] %in% data0[data0[, "SNP"] %in% eSNPlistRandom, "IDBLOCK"], "SNP"]
                }
                xSNPenrichStats <- table(assoc_xSNP, factor(xList%in%xSNPlistRandom, levels = c(FALSE, TRUE)))
                xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
            } else {
                xTmp <- NULL
            }
            return(list(eSNP = eTmp, xSNP = xTmp))
        })
        eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function(l){l$eSNP})))

        eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
        if (is.matrix(eSNPanyDup)) {
            if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
                eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
            } else {}
            if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
                eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[is.finite(eSNPsample[, 5]), 5])*1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
            } else {}
            eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
            Ze <- (eEnrichRatio-eMeanEnrichRatio)/((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/length(eSNPsample[, 5]))
            # Ze <- (eEnrichRatio-mean(eSNPsample[, 5]))/((sum((eSNPsample[, 5]-mean(eSNPsample[, 5]))^2))/length(eSNPsample[, 5]))
        } else {
            Ze <- as.numeric(NA)
        }
        
        if (isLD) {
            xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP}))) 
            xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
            if (is.matrix(xSNPanyDup)) {
                if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
                    xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
                } else {}
                if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
                    xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[is.finite(xSNPsample[, 5]), 5])*1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
                } else {}
                xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
                Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
                # Zx <- (xEnrichRatio-mean(xSNPsample[, 5]))/((sum((xSNPsample[, 5]-mean(xSNPsample[, 5]))^2))/length(xSNPsample[, 5]))
            } else {
                Zx <- as.numeric(NA)
            }
        } else {
            Zx <- as.numeric(NA)
        }
        iSample <- iSample + nSampleMin
        if ((iSample-(iSample%/%catStep*catStep)) == 0) {cat(iSample, ".. ", sep = "")}
        rm(resParallel)
    }

    eSNP <- newSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = pnorm(Ze, lower.tail = FALSE), Resampling = eSNPsample)
    if (isLD) {
        xSNP <- newSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = pnorm(Zx, lower.tail = FALSE), Resampling = xSNPsample)
    } else {
        xSNP <- newSNP()
    }
    object1@eSNP <- eSNP
    object1@xSNP <- xSNP
    return(object1)
}
