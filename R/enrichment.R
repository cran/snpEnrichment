######################################################################
########################## Class Enrichment ##########################
############################## Creation ##############################
######################################################################


setClass(
    Class = "Enrichment", 
    representation = representation(
        Loss = "data.frame", 
        Call = "list", 
        eSNP = "EnrichSNP", 
        xSNP = "EnrichSNP", 
        Chromosomes = "list"
    ), 
    prototype = prototype(
        Loss = data.frame(), 
        Call = list(readEnrichment = list(pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL, 
                                            distThresh = NULL, sigThresh = NULL, LD = NULL, ldDir = NULL, extendMethod = NULL, mc.cores = NULL), 
                        reSample = list(object = NULL, nSample = NULL, MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL)), 
        eSNP = enrichSNP(), 
        xSNP = enrichSNP(), 
        Chromosomes = eval(parse(text = paste0("list(", paste(paste0("Chrom", seq(22), " = chromosome()"), collapse = ", "), ")")))
    )
)


setMethod(f = "enrichment", signature = "ANY", definition = function(Loss, Call, eSNP, xSNP, Chromosomes){
    if (missing(Loss)) {
        Loss <- data.frame()
    } else {}
    if (missing(Call)) {
        Call <- list(readEnrichment = list(pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL, 
                                            distThresh = NULL, sigThresh = NULL, LD = NULL, ldDir = NULL, extendMethod = NULL, mc.cores = NULL), 
                        reSample = list(object = NULL, nSample = NULL, MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL))
    } else {}
    if (missing(Chromosomes)) {
        Chromosomes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq(22), " = chromosome()"), collapse = ", "), ")")))
    } else {}
    if (missing(eSNP)) {
        List <- eval(parse(text = paste0("c(", paste(paste0("Chromosomes$Chrom", seq(22), "@eSNP@List"), collapse=", "), ")")))
        eSNP <- enrichSNP(List = List)
    } else {}
    if (missing(xSNP)) {
        List <- eval(parse(text = paste0("c(", paste(paste0("Chromosomes$Chrom", seq(22), "@xSNP@List"), collapse=", "), ")")))
        xSNP <- enrichSNP(List = List)
    } else {}
    return(new("Enrichment", Loss = Loss, Call = Call, eSNP = eSNP, xSNP = xSNP, Chromosomes = Chromosomes))
})


setMethod(f = "is.enrichment", signature = "ANY", definition = function(object){
    if (length(object)>1) {
        return(sapply(object, is.enrichment))
    } else {
        if (class(object) == "Enrichment") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
})


setMethod(f = "summary", signature = "Enrichment", definition = function(object, extended = TRUE, complete = TRUE){
    if (missing(object)){
        stop('[Enrichment:summary] "object" is missing.', call. = FALSE)
        return(invisible())
    } else {}
    if (!is.logical(extended)){
        stop('[Enrichment:summary] "extended" must be logical.', call. = FALSE)
        return(invisible())
    } else {}
    if (!is.logical(complete)){
        stop('[Enrichment:summary] "complete" must be logical.', call. = FALSE)
        return(invisible())
    } else {}
    if (complete) {
        i <- "Complete"
    } else {
        i <- "Genome"
    }
    if (extended) {
        ExtendType <- c("eSNP", "xSNP")
    } else {
        ExtendType <- "eSNP"
    }
    switch(EXPR = i, 
        "Genome" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (type in ExtendType) {
                EnrichmentRatio <- eval(parse(text = paste0('object@', type, '@EnrichmentRatio')))
                Z <- eval(parse(text = paste0('object@', type, '@Z')))
                PValue <- eval(parse(text = paste0('object@', type, '@PValue')))
                Resampling <- eval(parse(text = paste0('nrow(object@', type, '@Resampling)')))
                Data <- object@Loss["Signal", ncol(object@Loss)]
                List <- eval(parse(text = paste0('length(object@', type, '@List)')))
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
        }, 
        "Complete" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (type in ExtendType) {
                EnrichmentRatio <- eval(parse(text = paste0('object@', type, '@EnrichmentRatio')))
                Z <- eval(parse(text = paste0('object@', type, '@Z')))
                PValue <- eval(parse(text = paste0('object@', type, '@PValue')))
                Resampling <- eval(parse(text = paste0('nrow(object@', type, '@Resampling)')))
                Data <- object@Loss["Signal", ncol(object@Loss)]
                List <- eval(parse(text = paste0('length(object@', type, '@List)')))
                resTmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                            if (length(Z)==0) {NA} else {Z}, 
                            if (length(PValue)==0) {NA} else {PValue}, 
                            if (length(Resampling)==0) {NA} else {Resampling}, 
                            if (length(Data)==0) {NA} else {Data}, 
                            if (length(List)==0) {NA} else {List}
                )
                for (iChr in seq(22)) {
                    EnrichmentRatio <- eval(parse(text = paste0('object@Chromosomes$Chrom', iChr, '@', type, '@EnrichmentRatio')))
                    Z <- eval(parse(text = paste0('object@Chromosomes$Chrom', iChr, '@', type, '@Z')))
                    PValue <- eval(parse(text = paste0('object@Chromosomes$Chrom', iChr, '@', type, '@PValue')))
                    Resampling <- eval(parse(text = paste0('nrow(object@Chromosomes$Chrom', iChr, '@', type, '@Resampling)')))
                    Data <- eval(parse(text = paste0('nrow(object@Chromosomes$Chrom', iChr, '@Data)')))
                    List <- eval(parse(text = paste0('length(object@Chromosomes$Chrom', iChr, '@', type, '@List)')))
                    tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                            if (length(Z)==0) {NA} else {Z}, 
                            if (length(PValue)==0) {NA} else {PValue}, 
                            if (length(Resampling)==0) {NA} else {Resampling}, 
                            if (length(Data)==0) {NA} else {Data}, 
                            if (length(List)==0) {NA} else {List}
                    )
                    resTmp <- rbind(resTmp, tmp)
                }
                rownames(resTmp) <- c("Genome", paste0("Chrom", seq(22)))
                colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                res[[type]] <- resTmp
            }
            if (extended) {
                return(res)
            } else {
                return(res[[ExtendType]])
            }
        }, 
        stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
    )
})


.Enrichment.show <- function(object){
    .showArgs <- function(args) {
        tmpArgs <- args[[names(args)]]
        types <- lapply(tmpArgs, class)
        res <- NULL
        for (iArg in names(tmpArgs)) {
            res <- c(res, paste0(iArg, 
                                ifelse(types[[iArg]]=="character", "=\"", "="), 
                                ifelse(types[[iArg]]=="NULL" | types[[iArg]]=="name", deparse(tmpArgs[[iArg]]), tmpArgs[[iArg]]), 
                                ifelse(types[[iArg]]=="character", "\"", "")))
        }
        paste0(names(args), '(', paste(res, collapse = paste0(", \n", paste(rep(" ", nchar(names(args))+5), collapse = ""))), ')\n')
    }
    cat(" ~ Loss :", paste0("(", paste(dim(object@Loss), collapse = "x"), ")"))
    if (nrow(object@Loss) == 0) {
        cat("\nNA")
    } else {
        cat("\n")
        print(object@Loss[seq(2), ], quote = FALSE)
        cat("   ..... .....\n")
    }       
    cat("\n ~ Call :\n   ", .showArgs(object@Call["readEnrichment"]))
    cat("\n   ", .showArgs(object@Call["reSample"]))
    cat("\n ~ eSNP :")
        .EnrichSNP.show(object@eSNP)
    cat("\n ~ xSNP :")
        .EnrichSNP.show(object@xSNP)
    nbVoid <- 0
    for (slot in names(object@Chromosomes)) {
        eval(parse(text = paste0("nbVoid <- nbVoid + as.numeric(identical(chromosome(), ", paste0("object@Chromosomes$", slot), "))")))
    }
    cat("\n ~ Chromosomes :", paste(nbVoid, length(names(object@Chromosomes)), sep = "/"), "empty Chromosomes")
    cat("\n")
    return(invisible())
}
setMethod(f = "show", signature = "Enrichment", definition = function(object){cat("    ~~~ Class:", class(object), "~~~ \n"); .Enrichment.show(object)})
setMethod(f = "print", signature = "Enrichment", definition = function(x, ...){cat("    ~~~ Class:", class(x), "~~~ \n"); .Enrichment.show(x)})


setMethod(f = "[", signature = "Enrichment", definition = function(x, i, j, drop){
    nbChr <- length(x@Chromosomes)
    if (missing(j)) {
        switch(EXPR = i, 
            "Loss" = {return(x@Loss)}, 
            "Data" = {
                resData <- mclapply2(seq(22), mc.cores = min(22, detectCores()), function(iChr) {
                    return(x@Chromosomes[[iChr]]@Data)
                })
                return(do.call("rbind", resData))
            }, 
            "LD" = {
                resLD <- mclapply2(seq(22), mc.cores = min(22, detectCores()), function(iChr) {
                    return(x@Chromosomes[[iChr]]@LD)
                })
                return(unlist(resLD))
            }, 
            "Call" = {return(x@Call)}, 
            "eSNP" = {return(x@eSNP)}, 
            "xSNP" = {return(x@xSNP)}, 
            "Table" = {
                res <- list(eSNP = NULL, xSNP = NULL)
                for (type in c("eSNP", "xSNP")) {
                    res[[type]] <- eval(parse(text = paste0('x@', type, '@Table')))
                }
                return(res)
            }, 
            "Chromosomes" = {return(x@Chromosomes)}, 
            "Stats" = {
                res <- list(eSNP = NULL, xSNP = NULL)
                for (type in c("eSNP", "xSNP")) {
                    EnrichmentRatio <- eval(parse(text = paste0('x@', type, '@EnrichmentRatio')))
                    Z <- eval(parse(text = paste0('x@', type, '@Z')))
                    PValue <- eval(parse(text = paste0('x@', type, '@PValue')))
                    Resampling <- eval(parse(text = paste0('nrow(x@', type, '@Resampling)')))
                    Data <- x@Loss["Signal", ncol(x@Loss)]
                    List <- eval(parse(text = paste0('length(x@', type, '@List)')))
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
                return(res)
            }, 
            stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
    } else {
        if (max(j)>nbChr) {
            if (j == "ALL") {
                switch(EXPR = i, 
                    "Loss" = {return(x@Loss)}, 
                    "Data" = {
                        resData <- mclapply2(seq(22), mc.cores = min(22, detectCores()), function(iChr) {
                            return(x@Chromosomes[[iChr]]@Data)
                        })
                        return(do.call("rbind", resData))
                    }, 
                    "LD" = {
                        resLD <- mclapply2(seq(22), mc.cores = min(22, detectCores()), function(iChr) {
                            return(x@Chromosomes[[iChr]]@LD)
                        })
                        return(unlist(resLD))
                    }, 
                    "Call" = {return(x@Call)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@List')))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@Table')))
                        }
                        return(res)
                    }, 
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@EnrichmentRatio')))
                        }
                        return(res)
                    }, 
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@Z')))
                        }
                        return(res)
                    }, 
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@PValue')))
                        }
                        return(res)
                    }, 
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('x@', type, '@Resampling')))
                        }
                        return(res)
                    }, 
                    "Chromosomes" = {return(x@Chromosomes)}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            EnrichmentRatio <- eval(parse(text = paste0('x@', type, '@EnrichmentRatio')))
                            Z <- eval(parse(text = paste0('x@', type, '@Z')))
                            PValue <- eval(parse(text = paste0('x@', type, '@PValue')))
                            Resampling <- eval(parse(text = paste0('nrow(x@', type, '@Resampling)')))
                            Data <- x@Loss["Signal", ncol(x@Loss)]
                            List <- eval(parse(text = paste0('length(x@', type, '@List)')))
                            resTmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                            )
                            for (iChr in seq(22)) {
                                EnrichmentRatio <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@EnrichmentRatio')))
                                Z <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@Z')))
                                PValue <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@PValue',)))
                                Resampling <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', iChr, '@', type, '@Resampling)')))
                                Data <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', iChr, '@Data)')))
                                List <- eval(parse(text = paste0('length(x@Chromosomes$Chrom', iChr, '@', type, '@List)')))
                                tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                                )
                                resTmp <- rbind(resTmp, tmp)
                            }
                            rownames(resTmp) <- c("Genome", paste0("Chrom", seq(22)))
                            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                            res[[type]] <- resTmp
                        }
                        return(res)
                    }, 
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
                )
            } else {
                stop('[Enrichment:get] "j" is out of limits.', call. = FALSE)
            }
        } else {
            if (length(j)>1) {
                switch(EXPR = i, 
                    "Signal" = {return(eval(parse(text = paste0('list(', paste(paste0("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]"), collapse = ", "), ')'))))}, 
                    "Data" = {
                        resData <- mclapply2(j, mc.cores = min(length(j), detectCores()), function(iChr) {
                            return(x@Chromosomes[[iChr]]@Data)
                        })
                        return(do.call("rbind", resData))
                    }, 
                    "LD" = {
                        resLD <- mclapply2(j, mc.cores = min(length(j), detectCores()), function(iChr) {
                            return(x@Chromosomes[[iChr]]@LD)
                        })
                        return(unlist(resLD))
                    }, 
                    "Call" = {return(x@Call)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('c(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@List"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('list(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@Table"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('c(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@EnrichmentRatio"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('c(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@Z"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('c(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@PValue"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0('list(', paste(paste0("x@Chromosomes$Chrom", j, "@", type,"@Resampling"), collapse = ", "), ')')))
                        }
                        return(res)
                    }, 
                    "Chromosomes" = {return(x@Chromosomes[j])}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            resTmp <- NULL
                            for (iChr in j) {
                                EnrichmentRatio <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@EnrichmentRatio')))
                                Z <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@Z')))
                                PValue <- eval(parse(text = paste0('x@Chromosomes$Chrom', iChr, '@', type, '@PValue')))
                                Resampling <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', iChr, '@', type, '@Resampling)')))
                                Data <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', iChr, '@Data)')))
                                List <- eval(parse(text = paste0('length(x@Chromosomes$Chrom', iChr, '@', type, '@List)')))
                                tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                                )
                                resTmp <- rbind(resTmp, tmp)
                            }
                            rownames(resTmp) <- c(paste0("Chrom", j))
                            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                            res[[type]] <- resTmp
                        }
                        return(res)
                    }, 
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
                )
            } else {
                switch(EXPR = i, 
                    "Signal" = {return(eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]"))))}, 
                    "Data" = {return(eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@Data"))))}, 
                    "LD" = {return(eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@LD"))))}, 
                    "Call" = {return(x@Call)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@List")))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@Table")))
                        }
                        return(res)
                    }, 
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@EnrichmentRatio")))
                        }
                        return(res)
                    }, 
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@Z")))
                        }
                        return(res)
                    }, 
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@PValue")))
                        }
                        return(res)
                    }, 
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste0("x@Chromosomes$Chrom", j, "@", type,"@Resampling")))
                        }
                        return(res)
                    }, 
                    "Chromosomes" = {return(x@Chromosomes[[j]])}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            EnrichmentRatio <- eval(parse(text = paste0('x@Chromosomes$Chrom', j, '@', type, '@EnrichmentRatio')))
                            Z <- eval(parse(text = paste0('x@Chromosomes$Chrom', j, '@', type, '@Z')))
                            PValue <- eval(parse(text = paste0('x@Chromosomes$Chrom', j, '@', type, '@PValue')))
                            Resampling <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', j, '@', type, '@Resampling)')))
                            Data <- eval(parse(text = paste0('nrow(x@Chromosomes$Chrom', j, '@Data)')))
                            List <- eval(parse(text = paste0('length(x@Chromosomes$Chrom', j, '@', type, '@List)')))
                            tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                    if (length(Z)==0) {NA} else {Z}, 
                                    if (length(PValue)==0) {NA} else {PValue}, 
                                    if (length(Resampling)==0) {NA} else {Resampling}, 
                                    if (length(Data)==0) {NA} else {Data}, 
                                    if (length(List)==0) {NA} else {List}
                            )
                            names(tmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                            res[[type]] <- tmp
                        }
                        return(res)
                    }, 
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
                )
            }
        }
    }
    return(invisible())
})


setMethod(f = "[<-", signature = "Enrichment", definition = function(x, i, j, value){
    nbChr <- length(x@Chromosomes)
    if (missing(j)) {
        switch(EXPR = i, 
            "Loss" = {x@Loss <- value}, 
            "Data" = {stop('[Enrichment:set] "Data" is not available for Set function.', call. = FALSE)}, 
            "LD" = {stop('[Enrichment:set] "LD" is not available for Set function.', call. = FALSE)}, 
            "Call" = {x@Call <- value}, 
            "eSNP" = {x@eSNP <- value}, 
            "xSNP" = {x@xSNP <- value}, 
            "List" = {stop('[Enrichment:set] "List" is not available for Set function.', call. = FALSE)}, 
            "Table" = {stop('[Enrichment:set] "Table" is not available for Set function.', call. = FALSE)}, 
            "EnrichmentRatio" = {stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function.', call. = FALSE)}, 
            "Z" = {stop('[Enrichment:set] "Z" is not available for Set function.', call. = FALSE)}, 
            "PValue" = {stop('[Enrichment:set] "PValue" is not available for Set function.', call. = FALSE)}, 
            "Resampling" = {stop('[Enrichment:set] "Resampling" is not available for Set function.', call. = FALSE)}, 
            "Chromosomes" = {x@Chromosomes <- value}, 
            "Stats" = {stop('[Enrichment:set] "Stats" is not available for Set function.', call. = FALSE)}, 
            stop('[Enrichment:set] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
        )
    } else {
        if (max(j)>nbChr) {
            stop('[Enrichment:set] "j" is out of limits.', call. = FALSE)
        } else {
            if (length(j)>1) {
                stop('[Enrichment:set] "j" must be atomic.', call. = FALSE)
            } else {
                switch(EXPR = i, 
                    "Loss" = {x@Loss <- value}, 
                    "Data" = {stop('[Enrichment:set] "Data" is not available for Set function.', call. = FALSE)}, 
                    "LD" = {stop('[Enrichment:set] "LD" is not available for Set function.', call. = FALSE)}, 
                    "Call" = {x@Call <- value}, 
                    "eSNP" = {x@Chromosomes[[j]]@eSNP <- value}, 
                    "xSNP" = {x@Chromosomes[[j]]@xSNP <- value}, 
                    "List" = {stop('[Enrichment:set] "List" is not available for Set function.', call. = FALSE)}, 
                    "Table" = {stop('[Enrichment:set] "Table" is not available for Set function.', call. = FALSE)}, 
                    "EnrichmentRatio" = {stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function.', call. = FALSE)}, 
                    "Z" = {stop('[Enrichment:set] "Z" is not available for Set function.', call. = FALSE)}, 
                    "PValue" = {stop('[Enrichment:set] "PValue" is not available for Set function.', call. = FALSE)}, 
                    "Resampling" = {stop('[Enrichment:set] "Resampling" is not available for Set function.', call. = FALSE)}, 
                    "Chromosomes" = {x@Chromosomes[[j]] <- value}, 
                    "Stats" = {stop('[Enrichment:set] "Stats" is not available for Set functions.', call. = FALSE)}, 
                    stop('[Enrichment:set] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
                )
            }
        }
    }
    validObject(x)
    return(x)
})


setMethod(f = "computeER", signature = "Enrichment", definition = function(object, sigThresh = 0.05, mc.cores = 1) {
    if (!missing(object)) {
        object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr){
            data <- chr@Data
            chrLD <- length(chr@LD)
            for (type in c("eSNP", "xSNP")) {
                if (!(chrLD == 0 & type == "xSNP")) {
                    snpEnrich <- table(factor(data[, "PVALUE"]<sigThresh, levels = c(FALSE, TRUE)), factor(data[, type], levels = c(0, 1)))
                    colnames(snpEnrich) <- c("otherSNP", type)
                    rownames(snpEnrich) <- eval(parse(text = paste0('c("P>=', sigThresh, '", "P<', sigThresh, '")')))
                    chr[type] <- enrichSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
                } else {}
            }
            return(chr)
        })
        for (type in c("eSNP", "xSNP")) {
            bigEnrichment <- matrix(0, nrow = 2, ncol = 2)
            for (jChr in seq(22)) {
                bigEnrichment <- eval(parse(text = paste0('bigEnrichment + object@Chromosomes$Chrom', jChr, '@', type, '@Table')))
            }
            object[type]@Table <- bigEnrichment
            object[type]@EnrichmentRatio <- .enrichmentRatio(bigEnrichment)
        }
        return(object)
    } else {
        stop('[Enrichment:computeER] "Enrichment" object is required.', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "reSample", signature = "Enrichment", definition = function(object, nSample = 100, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = FALSE) {
    if (!missing(object)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:reSample] nSample was increased to 10.', call. = FALSE)
        } else {}
        extendMethod <- object@Call$readEnrichment$extendMethod
        sigThresh <- object@Call$readEnrichment$sigThresh

        cat("########### Resample Enrichment ############\n")
        nSampleOld <- object@Call$reSample$nSample
        if (onlyGenome == FALSE) {
            listRes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq(22), " = NULL"), collapse = ", "), ")")))
            for (iChr in seq(22)) {
                cat("  Chromosome ", if (nchar(iChr) == 1) {paste0("0", iChr)} else {iChr}, ": ", sep = "")
                listRes[[iChr]] <- reSample(object = object@Chromosomes[[iChr]], nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
                cat("END\n")
            }
            object@Chromosomes <- listRes
            rm(listRes)
        } else {}

        cat("  Genome       : ")
        result <- .reSample(object = object, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        cat("END\n")
        rm(object)
        
        sysCall <- sys.call(sys.parent())
        params <- as.list(sysCall[-1])
        formal <- as.list(names(formals(as.character(sysCall))))
        names(formal) <- formal
        names(params) <- names(formal)[seq(length(params))]

        for (iArg in names(formal)) {
            if (identical(grep(iArg, names(params)), integer(0))) {
                if (grep(iArg, names(formal)) > length(params)) {
                    formal[[iArg]] <- eval(parse(text = formal[[iArg]]))
                } else {
                    argTmp <- params[[grep(iArg, names(formal))]]
                    if (is.character(argTmp) | isS4(eval(params[[grep(iArg, names(formal))]]))) {
                        formal[[iArg]] <- params[[grep(iArg, names(formal))]]
                    } else {
                        formal[[iArg]] <- eval(parse(text = params[[grep(iArg, names(formal))]]))
                    }
                }
            } else {
                formal[[iArg]] <- params[[iArg]]
            }
        }
        if (is.numeric(nSampleOld)) {
            formal$nSample <- nSampleOld + formal$nSample
        } else {}
        result@Call$reSample <- formal[c("object", "nSample", "MAFpool", "mc.cores", "onlyGenome")]

        nameObject <- deparse(params[[1]])
        assign(nameObject, result, inherits = TRUE, envir = parent.frame(2))
        cat("######## Resample Enrichment Done ##########\n")
        return(invisible())
    } else {
        stop('[Enrichment:reSample] "Enrichment" object is required.', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "excludeSNP", signature = "Enrichment", definition = function(object, excludeFile, mc.cores = 1) {
    if (missing(excludeFile)){
        warning('[Enrichment:excludeSNP] argument "excludeFile" is missing' , call. = FALSE)
        return(invisible())
    } else {
        cat("########## Exclude SNP list Start ##########\n")
        if (all(class(try(close(file(excludeFile)), silent = TRUE))!="try-error")) {
            eSNPexclude <- read.delim(file = excludeFile, header = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, stringsAsFactors = FALSE, sep = "\t")
        } else {
            eSNPexclude <- excludeFile
        }
        if (class(eSNPexclude) %in% c("matrix", "data.frame")) {
            if (ncol(eSNPexclude)>1) {
                eSNPexclude <- eSNPexclude[, 1]
            } else {}
        } else {}
        eSNPexclude <- unlist(eSNPexclude, use.names = FALSE)
        extendMethod <- object@Call$readEnrichment$extendMethod
        resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
            chrObject <- eval(parse(text = paste0('object@Chromosomes$Chrom', iChr)))
            temp <- chrObject@Data
            if (extendMethod == "ld") {
                xSNPexclude <- intersect(temp[, "SNP"], unique(c(eSNPexclude, chrObject@LD[names(chrObject@LD) %in% eSNPexclude])))
            } else {
                if (extendMethod == "block") {
                    temp0 <- temp[, c("SNP", "IDBLOCK")]
                    xSNPexclude <- temp0[temp0[, "IDBLOCK"] %in% temp0[temp0[, "SNP"] %in% eSNPexclude, "IDBLOCK"], "SNP"]
                } else {}
            }
            temp[temp[, "SNP"]%in%xSNPexclude, "eSNP"] <- 0
            temp[temp[, "SNP"]%in%xSNPexclude, "xSNP"] <- 0
            res <- chromosome(Data = temp, LD = chrObject@LD)
            return(res)
        })
        names(resParallel) <- paste0("Chrom", seq(22))
        result <- enrichment(Loss = object@Loss, 
                                Chromosomes = resParallel, 
                                Call = list(readEnrichment = object@Call$readEnrichment, 
                                                    reSample = list(object = NULL, nSample = NULL, MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL)))
        rm(resParallel)
        GC()

        result <- computeER(object = result, sigThresh = object@Call$readEnrichment$sigThresh, mc.cores = mc.cores)
        cat("########### Update SNP list END ############\n")
        for (type in c("eSNP", "xSNP")) {
            cat("   ", length(setdiff(object[type]@List, result[type]@List)), " SNPs are removed from", type, "list.\n")
        }
        result@Loss <- cbind(result@Loss, 
                                exclude = c(result@Loss["Signal", "CIS"], 
                                            length(result["List", seq(22)][["eSNP"]]), 
                                            sapply(seq(22), function(iChr){length(result["List", iChr][["eSNP"]])})))
        
        cat("########### Exclude SNP list END ###########\n")
        return(result)
    }
})


setMethod(f = "reset", signature = "Enrichment", definition = function(object, i) {
    switch(EXPR = i, 
        "Loss" = {
            object@Loss <- data.frame()
        }, 
        "Data" = {
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Data")
        }, 
        "LD" = {
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "LD")
        }, 
        "Call" = {
            object@Call <- list(readEnrichment = list(pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL, 
                                                        distThresh = NULL, sigThresh = NULL, LD = NULL, ldDir = NULL, extendMethod = NULL, mc.cores = NULL), 
                                reSample = list(object = NULL, nSample = NULL, MAFpool = NULL, mc.cores = NULL, onlyGenome = NULL))
        }, 
        "eSNP" = {
            object@eSNP <- enrichSNP()
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "eSNP")
        }, 
        "xSNP" = {
            object@xSNP <- enrichSNP()
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "xSNP")
        }, 
        "List" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@List <- character()
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "List")
        }, 
        "Table" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Table <- matrix(0, ncol = 2, nrow = 2)
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Table")
        }, 
        "EnrichmentRatio" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@EnrichmentRatio <- numeric()
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "EnrichmentRatio")
        }, 
        "Z" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Z <- numeric()
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Z")
        }, 
        "PValue" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@PValue <- numeric()
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "PValue")
        }, 
        "Resampling" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Resampling <- matrix(0, ncol = 5, nrow = 0)
            }
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Resampling")
        }, 
        "Chromosomes" = {
            object@Chromosomes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq(22), " = chromosome()"), collapse = ", "), ")")))
        }, 
        stop('[Enrichment:reset] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
    )
    return(object)
})


setMethod(f = "compareEnrichment", signature = "ANY", definition = function(object.x, object.y, pattern = "Chrom", nSample = 100, mc.cores = 1, onlyGenome = FALSE) {
    if (!missing(object.x) & !missing(object.y)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:compareEnrichment] "nSample" was increased to 10.', call. = FALSE)
        } else {}
        if (class(object.x) == "Enrichment" & class(object.y) == "Enrichment") {
            if (nrow(object.x["Data"]) == 0) {
                stop('[Enrichment:compareEnrichment] "Enrichment" data is empty for "object.x".', call. = FALSE)
                return(invisible())
            } else {}
            if (nrow(object.y["Data"]) == 0) {
                stop('[Enrichment:compareEnrichment] "Enrichment" data is empty for "object.y".', call. = FALSE)
                return(invisible())
            } else {}
        } else {
            stop('[Enrichment:compareEnrichment] "Enrichment" object is required.', call. = FALSE)
            return(invisible())
        }
        
        sigThresh.x <- object.x@Call$readEnrichment$sigThresh
        sigThresh.y <- object.y@Call$readEnrichment$sigThresh
        if (!identical(sigThresh.x, sigThresh.y)) {
            warning(paste0('[Enrichment:compareEnrichment] "sigThresh" differs from "object.x" to "object.y".\n         "object.x" parameter is taken: ', deparse(sigThresh.x)), call. = FALSE)
        } else {}
        sigThresh <- sigThresh.x
        
        MAFpool.x <- object.x@Call$reSample$MAFpool
        MAFpool.y <- object.y@Call$reSample$MAFpool
        if (!identical(MAFpool.x, MAFpool.y)) {
            warning(paste0('[Enrichment:compareEnrichment] "MAFpool" differs from "object.x" to "object.y".\n         "object.x" parameter is taken: ', deparse(MAFpool.x)), call. = FALSE)
        } else {}
        MAFpool <- eval(MAFpool.x)
        if (missing(MAFpool) | is.null(MAFpool)) {
            MAFpool <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5)
        } else{}
        
        extendMethod.x <- object.x@Call$readEnrichment$extendMethod
        extendMethod.y <- object.y@Call$readEnrichment$extendMethod
        if (!identical(extendMethod.x, extendMethod.y)) {
            warning(paste0('[Enrichment:compareEnrichment] "extendMethod" differs from "object.x" to "object.y".\n         "object.x" parameter is taken: ', deparse(extendMethod.x)), call. = FALSE)
        } else {}
        extendMethod <- extendMethod.x

        l1 <- object.x@eSNP@List
        l2 <- object.y@eSNP@List
        if (identical(l1, l2)) {
            stop('[Enrichment:compareEnrichment] Both lists are identical.', call. = FALSE)
            return(invisible())
        } else {}
        if (is.null(object.x@Call$reSample$nSample)) {
            reSample(object = object.x, nSample = nSample, MAFpool = MAFpool, mc.cores = mc.cores, onlyGenome = onlyGenome)
        } else {}
        if (is.null(object.y@Call$reSample$nSample)) {
            reSample(object = object.y, nSample = nSample, MAFpool = MAFpool, mc.cores = mc.cores, onlyGenome = onlyGenome)
        } else {}
        
        # enrichObject1 <- object.x
        # enrichObject2 <- object.y
        object.x <- reset(object.x, "Resampling")
        object.y <- reset(object.y, "Resampling")
        
        if (length(l1)<length(l2)) {
            object1 <- object.x
            object2 <- object.y
            namesRes <- c("Enrichment_1", "Enrichment_2", "PVALUE_1", "PVALUE_2", "PVALUE", "nSAMPLE", "SNP_1", "SNP_2")
        } else {
            object1 <- object.y
            object2 <- object.x
            namesRes <- c("Enrichment_2", "Enrichment_1", "PVALUE_2", "PVALUE_1", "PVALUE", "nSAMPLE", "SNP_2", "SNP_1")
        }
        rm(object.x, object.y)
        object2 <- reset(object2, "Data")
        object2 <- reset(object2, "LD")

        cat("############ Comparison Start ##############\n")
        if (onlyGenome == FALSE) {
            listRes <- eval(parse(text = paste0("list(", paste(paste0("Chrom", seq(22), " = NULL"), collapse = ", "), ")")))
            for (iChr in seq(22)) {
                # cat("  **** Chromosome", if (nchar(iChr) == 1) {paste0("0", iChr)} else {iChr}, "Comparison Start ****\n")
                cat("  Chromosome ", if (nchar(iChr) == 1) {paste0("0", iChr)} else {iChr}, ": ", sep = "")
                listRes[[iChr]] <- .compareEnrich(object1 = object1@Chromosomes[[iChr]], object2 = object2@Chromosomes[[iChr]], nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
                if (identical(sort(object1@Chromosomes[[iChr]]@eSNP@List), sort(object2@Chromosomes[[iChr]]@eSNP@List))) {
                    listRes[[iChr]] <- reset(listRes[[iChr]], "Z")
                    listRes[[iChr]] <- reset(listRes[[iChr]], "Resampling")
                    listRes[[iChr]]@eSNP@PValue <- as.numeric(NA)
                    listRes[[iChr]]@xSNP@PValue <- as.numeric(NA)
                } else {}
                cat("END\n")
            }
        } else {}

        cat("  Genome       : ")
        result <- .compareEnrich(object1 = object1, object2 = object2, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        if (onlyGenome == FALSE) {
            result@Chromosomes <- listRes
        } else {}
        cat("END\n")

        res <- list(eSNP = NULL, xSNP = NULL)
        summaryObj1 <- summary(object1)
        summaryObj2 <- summary(object2)
        summaryRes <- summary(result)
        res[["eSNP"]] <- cbind(summaryObj1[["eSNP"]][, 1], summaryObj2[["eSNP"]][, 1], summaryObj1[["eSNP"]][, 3], summaryObj2[["eSNP"]][, 3], summaryRes[["eSNP"]][, 3:4], summaryObj1[["eSNP"]][, 6], summaryObj2[["eSNP"]][, 6])
        res[["xSNP"]] <- cbind(summaryObj1[["xSNP"]][, 1], summaryObj2[["xSNP"]][, 1], summaryObj1[["xSNP"]][, 3], summaryObj2[["xSNP"]][, 3], summaryRes[["xSNP"]][, 3:4], summaryObj1[["xSNP"]][, 6], summaryObj2[["xSNP"]][, 6])
        colnames(res[["eSNP"]]) <- namesRes
        colnames(res[["xSNP"]]) <- namesRes
        
        cat("############# Comparison End ###############\n")
        return(res)
        # return(list(summary = res, object1 = enrichObject1, object2 = enrichObject2))
    } else {
        stop('[Enrichment:compareEnrichment] "Enrichment" object is required.', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "plot", signature = "Enrichment", definition = function(x, onlyGenome = FALSE, chr = NULL, types = c("eSNP", "xSNP"), ...){
    if (is.null(unlist(x@Call$reSample))) {
        stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
        return(invisible())
    } else {}
    if (!identical(x@Call$reSample$onlyGenome, onlyGenome) & x@Call$reSample$onlyGenome) {
        onlyGenome <- x@Call$reSample$onlyGenome
    } else {}
    if (onlyGenome) {
        for (type in types) {
            if (nrow(x[type]@Resampling)==0) {
                stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
                return(invisible())
            } else {}
        }
        if (!is.null(chr) & !missing(chr)) {
            if (nrow(x@Chromosomes[[chr]][type]@Resampling)==0) {
                stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
                return(invisible())
            } else {}
        } else {}
    } else {
        for (type in types) {
            if (any(c(Genome = nrow(x[type]@Resampling)==0, unlist(lapply(x@Chromosomes, function(x) {nrow(x[type]@Resampling)==0}))))) {
                stop('[Enrichment:plot] "reSample" have to be run before "plot".', call. = FALSE)
                return(invisible())
            } else {}
        }
    }

    if (missing(chr) | is.null(chr)) {
        matrixER <- list(eSNP = NULL, xSNP = NULL)
        for (type in types) {
            ER <- x[type]@EnrichmentRatio
            resample <- x[type]@Resampling[, 5]
            Z <- NULL
            ERs <- NULL
            size <- length(resample)
            if (size >= 1000) {
                interv <- unique(c(seq(from = 3, to = size, by = floor(size/1000)), size))
            } else {
                interv <- unique(c(seq(from = 3, to = size, by = 1), size))
            }
            for (k in interv) {
                Z <- c(Z, (ER-mean(resample[1:k], na.rm = TRUE))/var(resample[1:k], na.rm = TRUE))
                ERs <- resample
            }
            matrixER[[type]] <- Z
        }
        if (!onlyGenome) {
            for (j in seq(22)) {
                for (type in types) {
                    ER <- x@Chromosomes[[j]][type]@EnrichmentRatio
                    resample <- x@Chromosomes[[j]][type]@Resampling[, 5]
                    Z <- NULL
                    ERs <- NULL
                    size <- length(resample)
                    if (size >= 1000) {
                        interv <- unique(c(seq(from = 3, to = size, by = floor(size/1000)), size))
                    } else {
                        interv <- unique(c(seq(from = 3, to = size, by = 1), size))
                    }
                    for (k in interv) {
                        Z <- c(Z, (ER-mean(resample[1:k], na.rm = TRUE))/var(resample[1:k], na.rm = TRUE))
                        ERs <- resample
                    }
                    matrixER[[type]] <- cbind(matrixER[[type]], Z)
                    matrixER[[type]] <- replace(matrixER[[type]], grep(TRUE, is.infinite(matrixER[[type]])), NA)
                }
            }
        
            colors <- rainbow(22)
            par(mfrow = c(1, length(types)))
            for (type in types) {
                matrixER[[type]] <- apply(matrixER[[type]], 2, scale)
                plot(x = interv, y = matrixER[[type]][, 1], ylab = "Z (scale and center)", xlab = type, type = "l", ylim = range(na.exclude(matrixER[[type]])), ...)
                res <- sapply(seq(ncol(matrixER[[type]][, -1])), function(iER) {
                    lines(x = interv, y = matrixER[[type]][, iER+1], type = "l", ylim = range(na.exclude(matrixER[[type]])), col = colors[iER+1])
                })
                colnames(matrixER[[type]]) <- c("Genome", paste0("Chrom", seq(22)))
            }
        } else {
            colors <- rainbow(22)
            par(mfrow = c(1, length(types)))
            for (type in types) {
                matrixER[[type]] <- scale(matrixER[[type]])
                plot(x = interv, y = matrixER[[type]], ylab = "Z (scale and center)", xlab = type, type = "l", ylim = range(na.exclude(matrixER[[type]])), ...)
            }
        }
        return(invisible())
    } else {
        matrixER <- list(eSNP = NULL, xSNP = NULL)
        colors <- rainbow(length(chr))
        par(mfrow = c(1, length(types)))
        for (j in chr) {
            for (type in types) {
                ER <- x@Chromosomes[[j]][type]@EnrichmentRatio
                resample <- x@Chromosomes[[j]][type]@Resampling[, 5]
                Z <- NULL
                ERs <- NULL
                size <- length(resample)
                if (size >= 1000) {
                    interv <- unique(c(seq(from = 3, to = size, by = floor(size/1000)), size))
                } else {
                    interv <- unique(c(seq(from = 3, to = size, by = 1), size))
                }
                for (k in interv) {
                    Z <- c(Z, (ER-mean(resample[1:k], na.rm = TRUE))/var(resample[1:k], na.rm = TRUE))
                    ERs <- resample
                }
                matrixER[[type]] <- cbind(matrixER[[type]], Z)
                matrixER[[type]] <- replace(matrixER[[type]], grep(TRUE, is.infinite(matrixER[[type]])), NA)
            }
        }
        for (type in types) {
            matrixER[[type]] <- apply(matrixER[[type]], 2, scale)
            plot(x = interv, y = matrixER[[type]][, 1], ylab = "Z (scale and center)", xlab = type, type = "l", ylim = range(na.exclude(matrixER[[type]])), col = colors[1], ...)
            if (length(chr)>2) {
                res <- sapply(seq(ncol(matrixER[[type]][, -1])), function(iER) {
                    lines(x = interv, y = matrixER[[type]][, iER+1], type = "l", ylim = range(na.exclude(matrixER[[type]])), col = colors[iER+1])
                })
            } else {
                if (length(chr) == 2) {
                     lines(x = interv, y = matrixER[[type]][, 2], type = "l", ylim = range(na.exclude(matrixER[[type]])), col = colors[2])
                } else {}
            }
            colnames(matrixER[[type]]) <- paste0("Chrom", chr)
        }
        return(invisible())
    }
})
