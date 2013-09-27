setGeneric(name = "enrichSNP", def = function(List, Table, EnrichmentRatio, Z, PValue, Resampling){standardGeneric("enrichSNP")})
setGeneric(name = "chromosome", def = function(Data, LD, eSNP, xSNP){standardGeneric("chromosome")})
setGeneric(name = "enrichment", def = function(Loss, Call, eSNP, xSNP, Chromosomes){standardGeneric("enrichment")})
setGeneric(name = "doLDblock", def = function(object, mc.cores = 1){standardGeneric("doLDblock")})
setGeneric(name = "excludeSNP", def = function(object, excludeFile, mc.cores = 1){standardGeneric("excludeSNP")})
setGeneric(name = "computeER", def = function(object, sigThresh = 0.05, mc.cores = 1){standardGeneric("computeER")})
setGeneric(name = "reSample", def = function(object, nSample = 100, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE, ...){standardGeneric("reSample")})
setGeneric(name = "compareEnrichment", def = function(object.x, object.y, pattern = "Chrom", nSample = 100, mc.cores = 1, onlyGenome = FALSE){standardGeneric("compareEnrichment")})
setGeneric(name = "is.enrichment", def = function(object){standardGeneric("is.enrichment")})
setGeneric(name = "is.chromosome", def = function(object){standardGeneric("is.chromosome")})
setGeneric(name = "is.EnrichSNP", def = function(object){standardGeneric("is.EnrichSNP")})
setGeneric(name = "reset", def = function(object, i){standardGeneric("reset")})
# setGeneric(name = "summary", def = function(object, extended = TRUE, complete = TRUE){standardGeneric("summary")})

setMethod(f = "reSample", signature = "ANY", definition = function(object, nSample, sigThresh, MAFpool, extendMethod, mc.cores) {
    if (!(is.enrichment(object) & is.chromosome(object))){
        stop('[Method:reSample] not available for "', class(object), '" object', call. = FALSE)
        return(invisible())
    } else {}
})
setMethod(f = "excludeSNP", signature = "ANY", definition = function(object, excludeFile, mc.cores = 1) {
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
    # END <- unlist(regmatches(path, regexec(".$", path)))
    # START <- unlist(regmatches(path, regexec("^.", path)))
    # if (START != "/") {
        # path <- paste0("/", path)
    # } else {}
    # if (END != "/") {
        # path <- paste0(path, "/")
    # } else {}
    path <- gsub("/*$", "/", path)
    return(path)
}


maxCores <- function(mc.cores = 1) {
    if (Sys.info()[["sysname"]] != "Windows") {
        mc.cores.old <- mc.cores
        sysMemFree <- system("egrep 'MemFree' /proc/meminfo", intern = TRUE)
        sysMemAvailable <- 0.90*as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemFree))
        sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"))[8])
        mc.cores <- max(min(as.integer(mc.cores), floor(sysMemAvailable/sysProc)), 1)
        if (mc.cores != mc.cores.old) {
            warning(paste0('[mclapply2] To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".'), call. = FALSE)
        } else {}
    } else {}
    return(mc.cores)
}


mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, 
                      mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
    return(mclapply(X = X, FUN = FUN, ..., 
                    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent, 
                    mc.cores = maxCores(mc.cores), mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive))
}


.writeSignal <- function(pattern, snpInfoDir, signalFile) {
    tmpDir <- gsub("\\\\", "/", tempdir())
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

    chrom.bim <- read.delim(file = paste0(snpInfoDir, pattern, ".bim"), header = FALSE, stringsAsFactors = FALSE, 
                            colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""), 
                            check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "", "", "", ""))
    signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
    eval(parse(text = paste0('write.table(signal, file = "', tmpDir, '/snpEnrichment/', pattern, 
                                '.signal", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")')))
    return(invisible())
}


.readSignal <- function(pattern) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    signal <- read.delim(file = paste0(tmpDir, "/snpEnrichment/", pattern, ".signal"), header = TRUE, stringsAsFactors = FALSE, 
                            colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""), 
                            check.names = FALSE, strip.white = TRUE, col.names =  c("", "SNP", "PVALUE"))
    return(signal)
}


.writeFreq <- function(pattern, snpInfoDir) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    IN <- paste0(snpInfoDir, pattern)
    OUT <- paste0(tmpDir, "/snpEnrichment/", pattern)
    plinkData <- read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(pattern)[, "SNP"])
    plinkFreq <- col.summary(plinkData$genotypes)
    plinkFreq <- cbind(snp.name = rownames(plinkFreq), MAF = plinkFreq[, "MAF"])
    plinkRes <- merge(plinkData$map, plinkFreq, by = "snp.name")
    plinkRes <- plinkRes[, c("chromosome", "snp.name", "position", "MAF")]
    plinkRes[, "MAF"] <- as.numeric(as.character(plinkRes[, "MAF"]))
    colnames(plinkRes) <- c("CHR", "SNP", "POS", "MAF")
    write.table(plinkRes, paste0(OUT, ".all"), row.names = FALSE, sep = "\t")
    return(invisible())
}


writeLD <- function(pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, mc.cores = 1) {
    if (missing(pattern) | missing(snpInfoDir) | missing(signalFile) | missing(ldThresh)) {
        stop('[Enrichment:writeLD] argument(s) missing', call. = FALSE)
        return(invisible())
    } else {}
    tmpDir <- gsub("\\\\", "/", tempdir())
    dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
    snpInfoDir <- .checkFilePath(snpInfoDir)
    if (missing(ldDir) | is.null(ldDir)) {
        ldDir <- paste0(tmpDir, "/snpEnrichment/")
    } else {
        ldDir <- .checkFilePath(ldDir)
    }
    FILES <- list.files(snpInfoDir)
    cat("Compute LD for chromosome:\n  ")
    resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
        newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, ".*.bim"), FILES, value = TRUE), ".bim"))[1]
        isThereSignals <- grep(".signal", list.files(paste0(tmpDir, "/snpEnrichment/"), full.names = TRUE), value = TRUE)
        if (length(isThereSignals) != 22) {
            .writeSignal(pattern = newPattern, snpInfoDir, signalFile)
        } else {}
        IN <- paste0(snpInfoDir, newPattern)
        OUT <- paste0(tmpDir, "/snpEnrichment/", newPattern)
        plinkData <- read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(newPattern)[, "SNP"])
        ldData <- ld(x = plinkData$genotypes, depth = min(ncol(plinkData$genotypes)-1, 1000), stats = "R.squared")
        ldData <- as.matrix(ldData)
        ldData <- replace(ldData, grep(TRUE, is.na(ldData)), 0)
        ldData <- apply(ldData, 2, function(li) {which(li>0.8)})
        resLD <- matrix(unlist(strsplit(names(unlist(ldData)), "\\.")), ncol = 2, byrow = TRUE)
        write.table(resLD, file = paste0(ldDir, newPattern, ".ld"), sep = "", row.names = FALSE, col.names = FALSE)
        cat(iChr, " ", sep = "")
    })
    cat("\n")
    return(invisible())
}


initFiles <- function(pattern = "Chrom", snpInfoDir, snpListDir, signalFile, mc.cores = 1) {
    if (missing(snpInfoDir) | missing(signalFile)) {
        stop('[Enrichment:initFiles] argument(s) missing', call. = FALSE)
        return(invisible())
    } else {}
    snpInfoDir <- .checkFilePath(snpInfoDir)
    snpListDir <- .checkFilePath(snpListDir)
    FILES <- list.files(snpInfoDir)
    tmpDir <- gsub("\\\\", "/", tempdir())
    dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
    cat("All files are ready for chromosome:\n  ")
    resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
        newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, ".*.bim"), FILES, value = TRUE), ".bim"))[1]
        err <- try(.writeSignal(pattern = newPattern, snpInfoDir = snpInfoDir, signalFile = signalFile), silent = TRUE)
        warn <- warnings()
        if (class(err)=="try-error") {
            return(warn)
        } else {}
        err <- try(.writeFreq(pattern = newPattern, snpInfoDir = snpInfoDir), silent = TRUE)
        warn <- warnings()
        if (class(err)=="try-error") {
            return(warn)
        } else {}
        cat(iChr, " ", sep = "")
        return(invisible())
    })
    warn <- resParallel[[1]]
    if (!is.null(warn)) {
        stop(paste0("[Enrichment:initFiles] ", gsub("open file \'[^\']*'", paste0("create files in \'", snpInfoDir, "\'"), names(warn))), call. = FALSE)
    } else {}
    cat("\n")
    return(invisible())
}


.readSNP <- function(pattern, snpListDir) {
    snpListFile <- grep(pattern, list.files(snpListDir, full.names = TRUE), value = TRUE)[1]
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
    tmpDir <- gsub("\\\\", "/", tempdir())
    freq <- read.delim(file = paste0(tmpDir, "/snpEnrichment/", pattern, ".all"), header = TRUE, stringsAsFactors = FALSE, 
                        colClasses = c("numeric", "character", "numeric", "numeric"), na.string = c("NA", ""), 
                        check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "POS", "MAF"))
    return(freq)
}


.readLD <- function(pattern, snpInfoDir, ldDir) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    if (missing(ldDir) | is.null(ldDir)) {
        IN <- paste0(tmpDir, "/snpEnrichment/", pattern, ".ld")
    } else {
        IN <- paste0(.checkFilePath(ldDir), pattern, ".ld")
    }
    ldData <- read.delim(file = IN, header = FALSE, stringsAsFactors = FALSE, colClasses = c("character", "character"), 
                            na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = "")
    return(ldData)
}


.readFiles <- function(pattern, snpInfoDir, snpListDir, transcriptFile, distThresh) {
    newPattern <- unlist(strsplit(grep(paste0(pattern, ".*.bim"), list.files(snpInfoDir), value = TRUE), ".bim"))[1]
    eSNP <- .readSNP(pattern = newPattern, snpListDir = snpListDir)
    signal <- .readSignal(pattern = newPattern)
    plinkData <- .readFreq(pattern = newPattern, snpInfoDir = snpInfoDir)
    
    signalPlink <- merge(signal, plinkData, by = "SNP")
    if (nrow(eSNP)!=0) {
        eSNPunique <- eSNP[!duplicated(eSNP[, "SNP"]), ]
        snpSignal <- merge(signalPlink, eSNPunique[, c("SNP", "eSNP")], by = "SNP", all.x = TRUE)
        snpSignal[, "eSNP"][is.na(snpSignal[, "eSNP"])] <- 0
    } else {
        snpSignal <- signalPlink
        snpSignal[, "eSNP"] <- 0
    }
    
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
        signalLoss <- c(NA, length(unique(data[, "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])), sum(data[, "CIS"]))
        data <- data[data[, "CIS"] == 1, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")]
    } else {
        temp <- unique(snpSignal[, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")])
        data <- transform(temp, SNP = as.character(temp$SNP), PVALUE = as.numeric(temp$PVALUE), CHR = as.numeric(temp$CHR), POS = as.numeric(temp$POS), MAF = as.numeric(temp$MAF), eSNP = as.numeric(temp$eSNP))
        signalLoss <- c(NA, length(unique(data[, "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])), length(unique(data[!is.na(data[, "PVALUE"]), "SNP"])))
    }
    data[, "xSNP"] <- 0
    if (any(duplicated(data[, "SNP"]))) {
        cat(data[, "SNP"][duplicated(data[, "SNP"])], "\n")
        stop('[Enrichment:readEnrichment] Duplicated SNPs in Signal', call. = FALSE)
        return(invisible())
    } else {
        rownames(data) <- data[, "SNP"]
    }
    snpLoss <- c(length(eSNP[, "SNP"]), 
                    length(unique(eSNP[, "SNP"])), 
                    length(unique(intersect(eSNP[, "SNP"], signal[!is.na(signal[, "PVALUE"]), "SNP"]))), 
                    length(data[data[, "eSNP"]==1 & !is.na(data[, "PVALUE"]), "SNP"]))
    return(list(data = data, snpLoss = snpLoss, signalLoss = signalLoss))
}


readEnrichment <- function(pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 1000, sigThresh = 0.05, LD = FALSE, ldDir = NULL, extendMethod = "ld", mc.cores = 1) {
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
    resParallel <- mclapply2(seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
        files <- .readFiles(pattern = paste0(pattern, iChr), snpInfoDir = snpInfoDir, snpListDir = snpListDir, 
                            transcriptFile = transcriptFile, distThresh = distThresh)
        if (LD) {
            newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, ".*.bim"), list.files(snpInfoDir), value = TRUE), ".bim"))[1]
            linkageData <- .readLD(pattern = newPattern, snpInfoDir = snpInfoDir, ldDir = ldDir)
            
            dataTmp <- files$data
            data <- dataTmp[order(dataTmp$POS), ]
            snpSignal <- data[, "SNP"]
            
            linkageData <- linkageData[linkageData[, 1]%in%snpSignal & linkageData[, 2]%in%snpSignal, ]
            linkageData <- linkageData[order(data[linkageData[, 1], "POS"]), ]
            ldData <- c(linkageData[, 1], linkageData[, 2])
            names(ldData) <- c(linkageData[, 2], linkageData[, 1])
            if (extendMethod == "ld") {
                eSNP <- data[data[, "eSNP"]==1, "SNP"]
                xSNP <- intersect(data[, "SNP"], unique(c(data[data[, "eSNP"]==1, "SNP"], ldData[names(ldData) %in% eSNP])))
                data[intersect(xSNP, snpSignal), "xSNP"] <- 1
                data <- data[!is.na(data$PVALUE), ]
                resChr <- chromosome(Data = data, LD = ldData)
            } else {
                filesTmp <- files[c("snpLoss", "signalLoss")]
                rm(files, dataTmp, linkageData, snpSignal)
                GC()
                if (extendMethod == "block") { # extended LD for SNP which are not signal!!!!
                    resTmp <- chromosome(Data = data, LD = ldData)
                    resTmp <- doLDblock(object = resTmp, mc.cores = mc.cores)
                    data <- resTmp@Data
                    rm(resTmp)
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
            resChr <- chromosome(Data = files$data)
        }
        cat(iChr, " ", sep = "")
        return(list(resChr, files$snpLoss, files$signalLoss))
    })
    names(resParallel) <- paste0("Chrom", seq(22))
    cat("\n")
    
    result = enrichment()
    result@Chromosomes <- lapply(resParallel, function(l) {return(l[[1]])})
    snpLoss <- t(as.data.frame(lapply(resParallel, function(l) {return(l[[2]])})))
    signalLoss <- apply(t(as.data.frame(lapply(resParallel, function(l) {return(l[[3]])}))), 2, sum)
    loss <- rbind(Signal = signalLoss, Genome = apply(snpLoss, 2, sum), snpLoss)
    colnames(loss) <- c("Rows", "Unique", "Intersect.Signal", "CIS")
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

    sysCall <- sys.call(sys.parent())
    argsSNP <- as.list(sysCall[-1])
    formal <- as.list(names(formals(as.character(sysCall))))
    names(formal) <- formal

    names(argsSNP)[grep("^$", names(argsSNP))] <- names(formal)[grep("^$", names(argsSNP))]
    argsSNP <- c(argsSNP, lapply(formal[!names(formal) %in% names(argsSNP)], as.name))[names(formal)]
    paramClass <- sapply(argsSNP, class)

    for (iArg in names(formal)) {
        paramPos <- grep(iArg, names(formal))
        argTmp <- argsSNP[[paramPos]]
        classTmp <- paramClass[[paramPos]]
        switch(EXPR = classTmp,
            "character" = {formal[[iArg]] <- argTmp},
            "logical" = {formal[[iArg]] <- argTmp},
            "numeric" = {formal[[iArg]] <- argTmp},
            "integer" = {formal[[iArg]] <- argTmp},
            "NULL" = {formal[[iArg]] <- "NULL"},
            "call" = {formal[[iArg]] <- eval(argTmp)},
            "name" = {
                resEval <- eval(argTmp)
                switch(EXPR = class(resEval),
                    "character" = {formal[[iArg]] <- resEval},
                    "logical" = {formal[[iArg]] <- resEval},
                    "numeric" = {formal[[iArg]] <- resEval},
                    "integer" = {formal[[iArg]] <- resEval},
                    "matrix" = {formal[[iArg]] <- argTmp},
                    "data.frame" = {formal[[iArg]] <- argTmp}
                )
            }
        )
    }
    
    tmpDir <- gsub("\\\\", "/", tempdir())
    unlink(paste0(tmpDir, "/snpEnrichment"), recursive = TRUE)
    if (is.null(ldDir)) {
        formal[["ldDir"]] <- paste0(tmpDir, "/snpEnrichment/")
    } else {}
    result@Call$readEnrichment <- formal

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
        if (nResampling == 0) {
            xSNPsample <- NULL
        } else {
            xSNPsample <- object@xSNP@Resampling
        }
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
    if (isLD) {
        xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP})))
        xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
        Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
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
        if (isLD) {
            xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function(l){l$xSNP})))
            xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
            Zx <- (xEnrichRatio-xMeanEnrichRatio)/((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/length(xSNPsample[, 5]))
        } else {
            Zx <- 1
        }
        iSample <- iSample + nSampleMin
        if ((iSample-(iSample%/%catStep*catStep)) == 0) {cat(iSample, ".. ", sep = "")}
        rm(resParallel)
    }
    # cat("\n")


    object@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = pnorm(Ze, lower.tail = FALSE), Resampling = eSNPsample)
    if (isLD) {
        object@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = pnorm(Zx, lower.tail = FALSE), Resampling = xSNPsample)
    } else {
        object@xSNP <- enrichSNP()
    }
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
    object1@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = pnorm(Ze, lower.tail = FALSE), Resampling = eSNPsample)
    if (isLD) {
        object1@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = pnorm(Zx, lower.tail = FALSE), Resampling = xSNPsample)
    } else {
        object1@xSNP <- enrichSNP()
    }
    if (is.enrichment(object1)) {
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Z")
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "PValue")
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Resampling")
    } else {}
    return(object1)
}
