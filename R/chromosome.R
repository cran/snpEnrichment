######################################################################
########################## Class Chromosome ##########################
############################## Creation ##############################
######################################################################


setClass(
    Class = "Chromosome", 
    representation = representation(
        Data = "data.frame",
        LD = "character",
        eSNP = "SNP", 
        xSNP = "SNP" 
    ), 
    prototype = prototype(
        Data = data.frame(),
        LD = character(),
        eSNP = newSNP(), 
        xSNP = newSNP() 
    )
)


setMethod(f = "chromosome", signature = "ANY", definition = function(Data, LD, eSNP, xSNP){
    if (missing(Data)) {
        Data <- data.frame()
        if (missing(eSNP)) {eSNP <- xSNP <- newSNP()} else {}
    } else {
        if (missing(eSNP)) {
            eSNP <- newSNP(List = Data[Data[, "eSNP"] == 1, "SNP"])
            xSNP <- newSNP(List = Data[Data[, "xSNP"] == 1, "SNP"])
        } else {}
    }
    if (missing(LD)) {LD <- character()} else {}
    return(new("Chromosome", Data = Data, LD = LD, eSNP = eSNP, xSNP = xSNP))
})


setMethod(f = "is.chromosome", signature = "ANY", definition = function(object){
    if (length(object)>1) {
        return(sapply(object, is.chromosome))
    } else {
        if (class(object) == "Chromosome") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
})


setMethod(f = "summary", signature = "Chromosome", definition = function(object, extended = TRUE){
    if (missing(object)){
        stop('[Chromosome:summary] "object" is missing', call. = FALSE)
        return(invisible())
    } else {}
    if (extended) {
        ExtendType <- c("eSNP", "xSNP")
    } else {
        ExtendType <- "eSNP"
    }
    res <- list(eSNP = NULL, xSNP = NULL)
    for (type in ExtendType) {
        EnrichmentRatio <- eval(parse(text = paste('object@', type,'@EnrichmentRatio', sep = "")))
        Z <- eval(parse(text = paste('object@', type,'@Z', sep = "")))
        PValue <- eval(parse(text = paste('object@', type,'@PValue', sep = "")))
        Resampling <- eval(parse(text = paste('nrow(object@', type,'@Resampling)', sep = "")))
        Data <- eval(parse(text = paste('nrow(object@Data)', sep = "")))
        List <- eval(parse(text = paste('length(object@', type,'@List)', sep = "")))
        resTmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                    if (length(Z)==0) {NA} else {Z}, 
                    if (length(PValue)==0) {NA} else {PValue}, 
                    if (length(Resampling)==0) {NA} else {Resampling}, 
                    if (length(Data)==0) {NA} else {Data}, 
                    if (length(List)==0) {NA} else {List}
        )
        names(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
        res[[type]] <- resTmp
    }
    if (extended) {
        return(res)
    } else {
        return(res[[ExtendType]])
    }
})


.Chromosome.show <- function(object){
    cat("  ~ Data :", paste("(", paste(dim(object@Data), collapse = "x"), ")", sep=""))
        nrowShow <- min(5 , nrow(object@Data))
        ncolShow <- min(10 , ncol(object@Data))
        if (nrow(object@Data) != 0) {
            cat("\n")
            print(object@Data[seq(nrowShow), seq(ncolShow)], quote = FALSE)
            cat("   ..... .....")
        } else {
            cat(" NA\n")
        }
    # cat("\n  ~ LD :", paste("(", length(object@LD), ")", sep = ""), if (length(object@LD) <= 5) {if (length(object@LD) == 0) {"NA"} else {object@LD}} else {paste(paste(object@LD[seq(5)], collapse = " "), "...")})
    cat("\n  ~ eSNP :")
        .SNP.show(object@eSNP)
    cat("\n  ~ xSNP :")
        .SNP.show(object@xSNP)
    cat("\n")
    return(invisible())
}
setMethod(f = "show", signature = "Chromosome", definition = function(object){cat("     ~~~ Class:", class(object), "~~~\n"); .Chromosome.show(object)})


setMethod(f = "[", signature = "Chromosome", definition = function(x, i, j, drop){
    switch(EXPR = i, 
        "Data" = {return(x@Data)}, 
        "LD" = {return(x@LD)}, 
        "eSNP" = {return(x@eSNP)}, 
        "xSNP" = {return(x@xSNP)}, 
        stop('[Chromosome:get] ', i, ' is not a "Chromosome" slot', call. = FALSE)
    )
    return(invisible())
})


setMethod(f = "[<-", signature = "Chromosome", definition = function(x, i, j, value){
    switch(EXPR = i, 
        "Data" = {x@Data <- value}, 
        "LD" = {x@LD <- value}, 
        "eSNP" = {x@eSNP <- value}, 
        "xSNP" = {x@xSNP <- value}, 
        stop('[Chromosome:get] ', i, ' is not a "Chromosome" slot', call. = FALSE)
    )
    validObject(x)
    return(x)
})


setMethod(f = "computeER", signature = "Chromosome", definition = function(object, sigThresh = 0.05, mc.cores = detectCores()) {
    if (!missing(object)) {
        object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr){
            data <- chr@Data
            chrLD <- length(chr@LD)
            for (type in c("eSNP", "xSNP")) {
                if (!(chrLD == 0 & type == "xSNP")) {
                    snpEnrich <- table(factor(chr@Data[, "PVALUE"]<sigThresh, levels = c(FALSE, TRUE)), factor(chr@Data[, type], levels = c(0, 1)))
                    colnames(snpEnrich) <- c("otherSNP", type)
                    rownames(snpEnrich) <- eval(parse(text=paste('c("P>=', sigThresh, '", "P<', sigThresh, '")', sep="")))
                    chr[type] <- newSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
                } else {}
            }
            return(chr)
        })
        return(object)
    } else {
        stop('[Chromosome:computeER] "Chromosome" object is required', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "doLDblock", signature = "Chromosome", definition = function(object, mc.cores = detectCores()) {
    if (!missing(object)) {
        data <- object@Data
        jChr <- data[, "CHR"][1]
        cat("    - do LD block for chromosome", if (nchar(jChr) == 1) {paste("0", jChr, sep = "")} else {jChr}, "...\n")

        nbCORES <- mc.cores
        nbCores <- max(1, round((nbCORES-22)/22))
        chrLD <- object@LD
        
        byBlock <- split(chrLD, names(chrLD))
        byBlock <- unique(mclapply2(byBlock, mc.cores = nbCores, function(i) {names(i) <- NULL; return(i)}))
        LDblock <- do.call("rbind", unique(mclapply2(seq(length(byBlock)), mc.cores = nbCores, function(jBlock) {
            range(data[byBlock[[jBlock]], "POS"])
        })))
        
        names(byBlock) <- NULL
        LDblock <- LDblock[order(LDblock[, 1]), ]
        rm(chrLD, byBlock)
        GC()

        blockLim <- NULL
        for (iBlock in seq(nrow(LDblock))) {
            if (iBlock == 1) {
                POS <- LDblock[iBlock, ]
                blockLim <- rbind(blockLim, POS)
                jBlock <- 1
            } else {
                POS <- LDblock[iBlock, ]
                if (max(blockLim[nrow(blockLim), ])<min(POS)) {
                    blockLim <- rbind(blockLim, POS)
                } else {
                    blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
                }
                if (iBlock == nrow(LDblock)) {
                    iBlock <- iBlock + 1
                    if (max(blockLim[nrow(blockLim), ])<min(POS)) {
                        blockLim <- rbind(blockLim, POS)
                    } else {
                        blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
                    }
                } else {}
            }
        }
        rm(LDblock)
        GC()
        blockLim <- cbind(blockLim, seq(nrow(blockLim)))
        colnames(blockLim) <- c("MIN", "MAX", "IDBLOCK")
        rownames(blockLim) <- seq(nrow(blockLim))
        blockLim <- cbind(blockLim, LENGTH = NA)
        resParallel <- mclapply2(seq(nrow(blockLim)), mc.cores = nbCores, function(li){
            blockLim[li, "LENGTH"] <- as.numeric(blockLim[li, "MAX"])-as.numeric(blockLim[li, "MIN"])
            return(blockLim[li, ])
        })
        blockLim <- do.call("rbind", resParallel)

        data <- data[order(data[, "POS"]), ]
        data[, c("MIN", "MAX", "IDBLOCK", "LENGTH", "MAFmedian")] <- NA
        tmpChr <- mclapply2(seq(nrow(blockLim)), mc.cores = nbCores, function(i){
            m <- blockLim[i, ]
            interv <- seq(from = which(data[, "POS"] == m["MIN"]), to = which(data[, "POS"] == m["MAX"]))
            data[interv, c("MIN", "MAX", "IDBLOCK", "LENGTH")] <- matrix(rep(m, length(interv)), nrow = length(interv), byrow = TRUE)
            data[data[, "IDBLOCK"]%in%m["IDBLOCK"], "MAFmedian"] <- median(data[data[, "IDBLOCK"]%in%m["IDBLOCK"], "MAF"])
            return(data[interv, ])
        })

        data <- do.call("rbind", tmpChr)
        rownames(data) <- data[, "SNP"]
        object@Data <- data
        return(object)
    } else {
        stop('[Chromosome:doLDblock] "Chromosome" Object is required', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "reSample", signature = "Chromosome", definition = function(object, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()) {
    if (!missing(object)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:reSample] nSample was increased to 100', call. = FALSE)
        } else {}
        result <- .reSample(object = object, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        return(result)
    } else {
        stop('[Enrichment:reSample] "Enrichment" object is required', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "reset", signature = "Chromosome", definition = function(object, i) {
    switch(EXPR = i, 
        "Data" = {object@Data <- data.frame()}, 
        "LD" = {object@LD <- character()}, 
        "eSNP" = {object@eSNP <- newSNP()}, 
        "xSNP" = {object@xSNP <- newSNP()},
        "List" = {for (type in c("eSNP", "xSNP")) {object[type]@List <- character()}}, 
        "Table" = {for (type in c("eSNP", "xSNP")) {object[type]@Table <- matrix(0, ncol = 2, nrow = 2)}}, 
        "EnrichmentRatio" = {for (type in c("eSNP", "xSNP")) {object[type]@EnrichmentRatio <- numeric()}}, 
        "Z" = {for (type in c("eSNP", "xSNP")) {object[type]@Z <- numeric()}}, 
        "PValue" = {for (type in c("eSNP", "xSNP")) {object[type]@PValue <- numeric()}}, 
        "Resampling" = {for (type in c("eSNP", "xSNP")) {object[type]@Resampling <- matrix(0, ncol = 0, nrow = 0)}}, 
        stop('[Enrichment:reset] ', i, ' is not a "Enrichment" slot', call. = FALSE)
    )
    return(object)
})
