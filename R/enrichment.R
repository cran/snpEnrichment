######################################################################
########################## Class Enrichment ##########################
############################## Creation ##############################
######################################################################


setClass(
    Class = "Enrichment", 
    representation = representation(
        Signal = "data.frame", 
        Parameters = "list", 
        eSNP = "SNP",
        xSNP = "SNP",
        Chromosomes = "list"
    ), 
    prototype = prototype(
        Signal = data.frame(), 
        Parameters = list(read = list(output = NULL, pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, snpInfoDir = NULL, distThresh = NULL, sigThresh = NULL, LD = NULL, extendMethod = NULL), 
                            reSample = list(object = NULL, nSample = NULL, sigThresh = NULL, MAFpool = NULL)), 
        eSNP = newSNP(),
        xSNP = newSNP(), 
        Chromosomes = eval(parse(text = paste("list(", paste(paste("Chrom", seq(22), " = chromosome()",  sep = ""), collapse = ", "), ")", sep = "")))
    )
)


setMethod(f = "enrichment", signature = "ANY", definition = function(Signal, Parameters, eSNP, xSNP, Chromosomes){
    if (missing(Signal)) {Signal <- data.frame()} else {}
    if (missing(Parameters)) {
        Parameters <- list(read = list(output = NULL, pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, 
                                        snpInfoDir = NULL, distThresh = NULL, sigThresh = NULL, LD = NULL, extendMethod = NULL), 
                            reSample = list(object = NULL, nSample = NULL, sigThresh = NULL, MAFpool = NULL))
    } else {}
    if (missing(Chromosomes)) {Chromosomes <- eval(parse(text = paste("list(", paste(paste("Chrom", seq(22), " = chromosome()",  sep = ""), collapse = ", "), ")", sep = "")))} else {}
    if (missing(eSNP)) {
        List <- eval(parse(text=paste("c(", paste(paste("Chromosomes$Chrom", seq(22), "@eSNP@List", sep = ""), collapse=", "), ")", sep="")))
        eSNP <- newSNP(List = List)
    } else {}
    if (missing(xSNP)) {
        List <- eval(parse(text=paste("c(", paste(paste("Chromosomes$Chrom", seq(22), "@xSNP@List", sep = ""), collapse=", "), ")", sep="")))
        xSNP <- newSNP(List = List)
    } else {}
    return(new("Enrichment", Signal = Signal, Parameters = Parameters, eSNP = eSNP, xSNP = xSNP, Chromosomes = Chromosomes))
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
        stop('[Enrichment:summary] "object" is missing', call. = FALSE)
        return(invisible())
    } else {}
    if (complete) {i <- "Complete"} else {i <- "Genome"}
    if (extended) {
        ExtendType <- c("eSNP", "xSNP")
    } else {
        ExtendType <- "eSNP"
    }
    switch(EXPR = i, 
        "Genome" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (type in ExtendType) {
                EnrichmentRatio <- eval(parse(text = paste('object@', type,'@EnrichmentRatio', sep = "")))
                Z <- eval(parse(text = paste('object@', type,'@Z', sep = "")))
                PValue <- eval(parse(text = paste('object@', type,'@PValue', sep = "")))
                Resampling <- eval(parse(text = paste('nrow(object@', type,'@Resampling)', sep = "")))
                Data <- eval(parse(text = paste('sum(object@Signal$IN)', sep = "")))
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
        }, 
        "Complete" = {
            res <- list(eSNP = NULL, xSNP = NULL)
            for (type in ExtendType) {
                EnrichmentRatio <- eval(parse(text = paste('object@', type,'@EnrichmentRatio', sep = "")))
                Z <- eval(parse(text = paste('object@', type,'@Z', sep = "")))
                PValue <- eval(parse(text = paste('object@', type,'@PValue', sep = "")))
                Resampling <- eval(parse(text = paste('nrow(object@', type,'@Resampling)', sep = "")))
                Data <- eval(parse(text = paste('sum(object@Signal$IN)', sep = "")))
                List <- eval(parse(text = paste('length(object@', type,'@List)', sep = "")))
                resTmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                            if (length(Z)==0) {NA} else {Z}, 
                            if (length(PValue)==0) {NA} else {PValue}, 
                            if (length(Resampling)==0) {NA} else {Resampling}, 
                            if (length(Data)==0) {NA} else {Data}, 
                            if (length(List)==0) {NA} else {List}
                )
                for (iChr in seq(22)) {
                    EnrichmentRatio <- eval(parse(text = paste('object@Chromosomes$Chrom', iChr, '@', type,'@EnrichmentRatio', sep = "")))
                    Z <- eval(parse(text = paste('object@Chromosomes$Chrom', iChr, '@', type,'@Z', sep = "")))
                    PValue <- eval(parse(text = paste('object@Chromosomes$Chrom', iChr, '@', type,'@PValue', sep = "")))
                    Resampling <- eval(parse(text = paste('nrow(object@Chromosomes$Chrom', iChr, '@', type,'@Resampling)', sep = "")))
                    Data <- eval(parse(text = paste('nrow(object@Chromosomes$Chrom', iChr, '@Data)', sep = "")))
                    List <- eval(parse(text = paste('length(object@Chromosomes$Chrom', iChr, '@', type,'@List)', sep = "")))
                    tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                            if (length(Z)==0) {NA} else {Z}, 
                            if (length(PValue)==0) {NA} else {PValue}, 
                            if (length(Resampling)==0) {NA} else {Resampling}, 
                            if (length(Data)==0) {NA} else {Data}, 
                            if (length(List)==0) {NA} else {List}
                    )
                    resTmp <- rbind(resTmp, tmp)
                }
                rownames(resTmp) <- c("Genome", paste("Chrom", seq(22), sep = ""))
                colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                res[[type]] <- resTmp
            }
            if (extended) {
                return(res)
            } else {
                return(res[[ExtendType]])
            }
        },
        stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot', call. = FALSE)
    )
})


.Enrichment.show <- function(object){
    cat(" ~ Signal :", paste("(", paste(dim(object@Signal), collapse = "x"), ")", sep = ""))
    cat("\n  * ", sum(object@Signal$IN == 0), " (", round((sum(object@Signal$IN == 0)/nrow(object@Signal))*100, digits = 2), " %) SNPs are removed", sep = "")
    if (!is.null(object@Parameters$read$sigThresh)) {
        cat("\n")
        distrib <- table(factor(object@Signal$PVALUE<object@Parameters$read$sigThresh, levels = c(FALSE, TRUE)), factor(object@Signal$IN, levels = c(0, 1)))
        rownames(distrib) <- c(paste("P>=", object@Parameters$read$sigThresh, sep = ""), paste("P<", object@Parameters$read$sigThresh, sep = ""))
        colnames(distrib) <- c("OUT", "IN")
        print(distrib, quote = FALSE)
    } else {}
    # nrowShow <- min(5 , nrow(object@Signal))
    # ncolShow <- min(5 , ncol(object@Signal))
    # if (nrow(object@Signal) == 0) {
        # cat("\nNA")
    # } else {
        # cat("\n")
        # print(object@Signal[seq(nrowShow), seq(ncolShow)], quote = FALSE)
        # cat("   ..... .....\n")
    # }       
    # cat("\n ~ Parameters :", .showArgs(object@Parameters["read"]))
    # cat("\n               ", .showArgs(object@Parameters["reSample"]))
    cat("\n ~ eSNP :")
        .SNP.show(object@eSNP)
    cat("\n ~ xSNP :")
        .SNP.show(object@xSNP)
    nbVoid <- 0
    for (slot in names(object@Chromosomes)) {
        eval(parse(text = paste("nbVoid <- nbVoid + as.numeric(identical(chromosome(), ", paste("object@Chromosomes$", slot, sep = ""), "))", sep = "")))
    }
    cat("\n ~ Chromosomes :", paste(nbVoid, length(names(object@Chromosomes)), sep = "/"), "empty Chromosomes")
    cat("\n")
    return(invisible())
}
setMethod(f = "show", signature = "Enrichment", definition = function(object){cat("    ~~~ Class:", class(object), "~~~ \n"); .Enrichment.show(object)})


setMethod(f = "[", signature = "Enrichment", definition = function(x, i, j, drop){
    nbChr <- length(x@Chromosomes)
    if (missing(j)) {
        switch(EXPR = i, 
            "Signal" = {return(x@Signal)}, 
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
            "Parameters" = {return(x@Parameters)}, 
            "eSNP" = {return(x@eSNP)}, 
            "xSNP" = {return(x@xSNP)}, 
            "Table" = {
                res <- list(eSNP = NULL, xSNP = NULL)
                for (type in c("eSNP", "xSNP")) {
                    res[[type]] <- eval(parse(text = paste('x@', type,'@Table', sep = "")))
                }
                return(res)
            },
            "Chromosomes" = {return(x@Chromosomes)}, 
            "Stats" = {
                res <- list(eSNP = NULL, xSNP = NULL)
                for (type in c("eSNP", "xSNP")) {
                    EnrichmentRatio <- eval(parse(text = paste('x@', type,'@EnrichmentRatio', sep = "")))
                    Z <- eval(parse(text = paste('x@', type,'@Z', sep = "")))
                    PValue <- eval(parse(text = paste('x@', type,'@PValue', sep = "")))
                    Resampling <- eval(parse(text = paste('nrow(x@', type,'@Resampling)', sep = "")))
                    Data <- eval(parse(text = paste('sum(x@Signal$IN)', sep = "")))
                    List <- eval(parse(text = paste('length(x@', type,'@List)', sep = "")))
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
            stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot', call. = FALSE)
        )
    } else {
        if (max(j)>nbChr) {
            if (j == "ALL") {
                switch(EXPR = i, 
                    "Signal" = {return(x@Signal)}, 
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
                    "Parameters" = {return(x@Parameters)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@List', sep = "")))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@Table', sep = "")))
                        }
                        return(res)
                    },
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@EnrichmentRatio', sep = "")))
                        }
                        return(res)
                    },
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@Z', sep = "")))
                        }
                        return(res)
                    },
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@PValue', sep = "")))
                        }
                        return(res)
                    },
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('x@', type,'@Resampling', sep = "")))
                        }
                        return(res)
                    },
                    "Chromosomes" = {return(x@Chromosomes)}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            EnrichmentRatio <- eval(parse(text = paste('x@', type,'@EnrichmentRatio', sep = "")))
                            Z <- eval(parse(text = paste('x@', type,'@Z', sep = "")))
                            PValue <- eval(parse(text = paste('x@', type,'@PValue', sep = "")))
                            Resampling <- eval(parse(text = paste('nrow(x@', type,'@Resampling)', sep = "")))
                            Data <- eval(parse(text = paste('sum(x@Signal$IN)', sep = "")))
                            List <- eval(parse(text = paste('length(x@', type,'@List)', sep = "")))
                            resTmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                            )
                            for (iChr in seq(22)) {
                                EnrichmentRatio <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@EnrichmentRatio', sep = "")))
                                Z <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@Z', sep = "")))
                                PValue <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@PValue', sep = "")))
                                Resampling <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', iChr, '@', type,'@Resampling)', sep = "")))
                                Data <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', iChr, '@Data)', sep = "")))
                                List <- eval(parse(text = paste('length(x@Chromosomes$Chrom', iChr, '@', type,'@List)', sep = "")))
                                tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                                )
                                resTmp <- rbind(resTmp, tmp)
                            }
                            rownames(resTmp) <- c("Genome", paste("Chrom", seq(22), sep = ""))
                            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                            res[[type]] <- resTmp
                        }
                        return(res)
                    }, 
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot', call. = FALSE)
                )
            } else {
                stop('[Enrichment:get] "j" is out of limits', call. = FALSE)
            }
        } else {
            if (length(j)>1) {
                switch(EXPR = i, 
                    "Signal" = {return(eval(parse(text = paste('list(', paste(paste("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]", sep = ""), collapse = ", "), ')', sep = ""))))}, 
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
                    "Parameters" = {return(x@Parameters)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('c(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@List", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('list(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@Table", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('c(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@EnrichmentRatio", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('c(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@Z", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('c(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@PValue", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste('list(', paste(paste("x@Chromosomes$Chrom", j, "@", type,"@Resampling", sep = ""), collapse = ", "), ')', sep = "")))
                        }
                        return(res)
                    }, 
                    "Chromosomes" = {return(x@Chromosomes[j])}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            resTmp <- NULL
                            for (iChr in j) {
                                EnrichmentRatio <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@EnrichmentRatio', sep = "")))
                                Z <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@Z', sep = "")))
                                PValue <- eval(parse(text = paste('x@Chromosomes$Chrom', iChr, '@', type,'@PValue', sep = "")))
                                Resampling <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', iChr, '@', type,'@Resampling)', sep = "")))
                                Data <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', iChr, '@Data)', sep = "")))
                                List <- eval(parse(text = paste('length(x@Chromosomes$Chrom', iChr, '@', type,'@List)', sep = "")))
                                tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
                                        if (length(Z)==0) {NA} else {Z}, 
                                        if (length(PValue)==0) {NA} else {PValue}, 
                                        if (length(Resampling)==0) {NA} else {Resampling}, 
                                        if (length(Data)==0) {NA} else {Data}, 
                                        if (length(List)==0) {NA} else {List}
                                )
                                resTmp <- rbind(resTmp, tmp)
                            }
                            rownames(resTmp) <- c(paste("Chrom", j, sep = ""))
                            colnames(resTmp) <- c("EnrichmentRatio", "Z", "PVALUE", "nbSample", "SNP", type)
                            res[[type]] <- resTmp
                        }
                        return(res)
                    }, 
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot', call. = FALSE)
                )
            } else {
                switch(EXPR = i, 
                    "Signal" = {return(eval(parse(text = paste("x@Chromosomes$Chrom", j, "@Data[, c('SNP', 'PVALUE')]", sep = ""))))}, 
                    "Data" = {return(eval(parse(text = paste("x@Chromosomes$Chrom", j, "@Data", sep = ""))))}, 
                    "LD" = {return(eval(parse(text = paste("x@Chromosomes$Chrom", j, "@LD", sep = ""))))}, 
                    "Parameters" = {return(x@Parameters)}, 
                    "List" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@List", sep = "")))
                        }
                        return(res)
                    }, 
                    "Table" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@Table", sep = "")))
                        }
                        return(res)
                    }, 
                    "EnrichmentRatio" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@EnrichmentRatio", sep = "")))
                        }
                        return(res)
                    }, 
                    "Z" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@Z", sep = "")))
                        }
                        return(res)
                    }, 
                    "PValue" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@PValue", sep = "")))
                        }
                        return(res)
                    }, 
                    "Resampling" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            res[[type]] <- eval(parse(text = paste("x@Chromosomes$Chrom", j, "@", type,"@Resampling", sep = "")))
                        }
                        return(res)
                    }, 
                    "Chromosomes" = {return(x@Chromosomes[[j]])}, 
                    "Stats" = {
                        res <- list(eSNP = NULL, xSNP = NULL)
                        for (type in c("eSNP", "xSNP")) {
                            EnrichmentRatio <- eval(parse(text = paste('x@Chromosomes$Chrom', j, '@', type,'@EnrichmentRatio', sep = "")))
                            Z <- eval(parse(text = paste('x@Chromosomes$Chrom', j, '@', type,'@Z', sep = "")))
                            PValue <- eval(parse(text = paste('x@Chromosomes$Chrom', j, '@', type,'@PValue', sep = "")))
                            Resampling <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', j, '@', type,'@Resampling)', sep = "")))
                            Data <- eval(parse(text = paste('nrow(x@Chromosomes$Chrom', j, '@Data)', sep = "")))
                            List <- eval(parse(text = paste('length(x@Chromosomes$Chrom', j, '@', type,'@List)', sep = "")))
                            tmp <- c(if (length(EnrichmentRatio)==0) {NA} else {EnrichmentRatio}, 
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
                    stop('[Enrichment:get] ', i, ' is not a "Enrichment" slot', call. = FALSE)
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
            "Signal" = {x@Signal <- value}, 
            "Data" = {stop('[Enrichment:set] "Data" is not available for Set function', call. = FALSE)}, 
            "LD" = {stop('[Enrichment:set] "LD" is not available for Set function', call. = FALSE)}, 
            "Parameters" = {x@Parameters <- value}, 
            "eSNP" = {x@eSNP <- value}, 
            "xSNP" = {x@xSNP <- value},
            "List" = {stop('[Enrichment:set] "List" is not available for Set function', call. = FALSE)}, 
            "Table" = {stop('[Enrichment:set] "Table" is not available for Set function', call. = FALSE)}, 
            "EnrichmentRatio" = {stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function', call. = FALSE)}, 
            "Z" = {stop('[Enrichment:set] "Z" is not available for Set function', call. = FALSE)}, 
            "PValue" = {stop('[Enrichment:set] "PValue" is not available for Set function', call. = FALSE)}, 
            "Resampling" = {stop('[Enrichment:set] "Resampling" is not available for Set function', call. = FALSE)}, 
            "Chromosomes" = {x@Chromosomes <- value}, 
            "Stats" = {stop('[Enrichment:set] "Stats" is not available for Set function', call. = FALSE)}, 
            stop('[Enrichment:set] ', i, ' is not a "Enrichment" slot', call. = FALSE)
        )
    } else {
        if (max(j)>nbChr) {
            stop('[Enrichment:set] "j" is out of limits', call. = FALSE)
        } else {
            if (length(j)>1) {
                stop('[Enrichment:set] "j" must be atomic', call. = FALSE)
            } else {
                switch(EXPR = i, 
                    "Signal" = {x@Signal <- value}, 
                    "Data" = {stop('[Enrichment:set] "Data" is not available for Set function', call. = FALSE)}, 
                    "LD" = {stop('[Enrichment:set] "LD" is not available for Set function', call. = FALSE)}, 
                    "Parameters" = {x@Parameters <- value}, 
                    "eSNP" = {x@Chromosomes[[j]]@eSNP <- value}, 
                    "xSNP" = {x@Chromosomes[[j]]@xSNP <- value},
                    "List" = {stop('[Enrichment:set] "List" is not available for Set function', call. = FALSE)}, 
                    "Table" = {stop('[Enrichment:set] "Table" is not available for Set function', call. = FALSE)}, 
                    "EnrichmentRatio" = {stop('[Enrichment:set] "EnrichmentRatio" is not available for Set function', call. = FALSE)}, 
                    "Z" = {stop('[Enrichment:set] "Z" is not available for Set function', call. = FALSE)}, 
                    "PValue" = {stop('[Enrichment:set] "PValue" is not available for Set function', call. = FALSE)}, 
                    "Resampling" = {stop('[Enrichment:set] "Resampling" is not available for Set function', call. = FALSE)}, 
                    "Chromosomes" = {x@Chromosomes[[j]] <- value}, 
                    "Stats" = {stop('[Enrichment:set] "Stats" is not available for Set functions', call. = FALSE)}, 
                    stop('[Enrichment:set] ', i, ' is not a "Enrichment" slot', call. = FALSE)
                )
            }
        }
    }
    validObject(x)
    return(x)
})


setMethod(f = "computeER", signature = "Enrichment", definition = function(object, sigThresh = 0.05, mc.cores = detectCores()) {
    if (!missing(object)) {
        object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function(chr){
            data <- chr@Data
            chrLD <- length(chr@LD)
            for (type in c("eSNP", "xSNP")) {
                if (!(chrLD == 0 & type == "xSNP")) {
                    snpEnrich <- table(factor(data[, "PVALUE"]<sigThresh, levels = c(FALSE, TRUE)), factor(data[, type], levels = c(0, 1)))
                    colnames(snpEnrich) <- c("otherSNP", type)
                    rownames(snpEnrich) <- eval(parse(text=paste('c("P>=', sigThresh, '", "P<', sigThresh, '")', sep="")))
                    chr[type] <- newSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
                } else {}
            }
            return(chr)
        })
        for (type in c("eSNP", "xSNP")) {
            bigEnrichment <- matrix(0, nrow = 2, ncol = 2)
            for (jChr in seq(22)) {
                bigEnrichment <- eval(parse(text = paste('bigEnrichment + object@Chromosomes$Chrom', jChr, '@', type, '@Table', sep = "")))
            }
            object[type]@Table <- bigEnrichment
            object[type]@EnrichmentRatio <- .enrichmentRatio(bigEnrichment)
        }
        return(object)
    } else {
        stop('[Enrichment:computeER] "Enrichment" object is required', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "reSample", signature = "Enrichment", definition = function(object, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()) {
    if (!missing(object)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:reSample] nSample was increased to 100', call. = FALSE)
        } else {}
        nSampleOld <- object@Parameters$reSample$nSample
        listRes <- eval(parse(text = paste("list(", paste(paste("Chrom", seq(22), " = NULL",  sep = ""), collapse = ", "), ")", sep = "")))
        for (iChr in seq(22)) {
            cat("  **** Chromosome", if (nchar(iChr) == 1) {paste("0", iChr, sep = "")} else {iChr}, "Computation Start ***\n")
            listRes[[iChr]] <- reSample(object = object@Chromosomes[[iChr]], nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
            cat("  ***** Chromosome", if (nchar(iChr) == 1) {paste("0", iChr, sep = "")} else {iChr}, "Computation End ****\n")
        }
        object@Chromosomes <- listRes
        rm(listRes)

        cat("  ****** Genome Computation Start  *******\n")
        result <- .reSample(object = object, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        cat("  ******* Genome Computation End  ********\n")
        rm(object)
        
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
                    if (is.character(argTmp)) {
                        formal[[iArg]] <- argsSNP[[grep(iArg, names(formal))]]
                    } else {
                        formal[[iArg]] <- eval(parse(text = argsSNP[[grep(iArg, names(formal))]]))
                    }
                }
            } else {
                formal[[iArg]] <- argsSNP[[iArg]]
            }
        }
        if (is.numeric(nSampleOld)) {
            formal$nSample <- nSampleOld + formal$nSample
        } else {}
        result@Parameters$reSample <- formal
    
        nameObject <- deparse(argsSNP[[1]])
        assign(nameObject, result, envir = parent.frame())
        return(invisible())
    } else {
        stop('[Enrichment:reSample] "Enrichment" object is required', call. = FALSE)
        return(invisible())
    }
})


setMethod(f = "excludeSNP", signature = "Enrichment", definition = function(object, excludeFile, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()) {
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
        
        resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function(iChr){
            chrObject <- eval(parse(text = paste('object@Chromosomes$Chrom', iChr, sep = "")))
            temp <- chrObject@Data
            temp[temp[, "SNP"]%in%eSNPexclude, "eSNP"] <- 0
            temp[temp[, "SNP"]%in%eSNPexclude, "xSNP"] <- 0
            res <- chromosome(Data = temp, LD = chrObject@LD)
            return(res)
        })
        names(resParallel) <- paste("Chrom", seq(22), sep = "")
        result <- enrichment(Signal = object@Signal, 
                                Chromosomes = resParallel, 
                                Parameters = list(read = object@Parameters$read, reSample = NULL))
        rm(resParallel)
        GC()

        result <- computeER(object = result, sigThresh = object@Parameters$read$sigThresh, mc.cores = mc.cores)
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:reSample] nSample was increased to 100', call. = FALSE)
        } else {}
        reSample(object = result, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        cat("########### Update SNP list END ############\n")
        cat("########### Exclude SNP list END ###########\n")
        return(result)
    }
})


setMethod(f = "reset", signature = "Enrichment", definition = function(object, i) {
    switch(EXPR = i, 
        "Signal" = {object@Signal <- data.frame()}, 
        "Data" = {object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Data")}, 
        "LD" = {object@Chromosomes <- lapply(object@Chromosomes, reset, i = "LD")}, 
        "Parameters" = {object@Parameters <- list(read = list(output = NULL, pattern = NULL, signalFile = NULL, transcriptFile = NULL, snpListDir = NULL, 
                                                            snpInfoDir = NULL, distThresh = NULL, sigThresh = NULL, LD = NULL, extendMethod = NULL), 
                                                    reSample = list(object = NULL, nSample = NULL, sigThresh = NULL, MAFpool = NULL))}, 
        "eSNP" = {object@eSNP <- newSNP()
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "eSNP")}, 
        "xSNP" = {object@xSNP <- newSNP()
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "xSNP")},
        "List" = {for (type in c("eSNP", "xSNP")) {object[type]@List <- character()}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "List")}, 
        "Table" = {for (type in c("eSNP", "xSNP")) {object[type]@Table <- matrix(0, ncol = 2, nrow = 2)}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Table")}, 
        "EnrichmentRatio" = {for (type in c("eSNP", "xSNP")) {object[type]@EnrichmentRatio <- numeric()}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "EnrichmentRatio")}, 
        "Z" = {for (type in c("eSNP", "xSNP")) {object[type]@Z <- numeric()}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Z")}, 
        "PValue" = {for (type in c("eSNP", "xSNP")) {object[type]@PValue <- numeric()}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "PValue")}, 
        "Resampling" = {for (type in c("eSNP", "xSNP")) {object[type]@Resampling <- matrix(0, ncol = 0, nrow = 0)}
            object@Chromosomes <- lapply(object@Chromosomes, reset, i = "Resampling")}, 
        "Chromosomes" = {object@Chromosomes <- eval(parse(text = paste("list(", paste(paste("Chrom", seq(22), " = chromosome()",  sep = ""), collapse = ", "), ")", sep = "")))}, 
        stop('[Enrichment:reset] ', i, ' is not a "Enrichment" slot', call. = FALSE)
    )
    return(object)
})


setMethod(f = "compareEnrichment", signature = "ANY", definition = function(list1, list2, pattern = "Chrom", signalFile, transcriptFile = FALSE, snpInfoDir, distThresh = 1000, LD = FALSE, ldThresh = 0.8, nSample = 100, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), extendMethod = "ld", mc.cores = detectCores()) {
    if (!missing(list1) & !missing(list2)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:compareEnrichment] nSample was increased to 100', call. = FALSE)
        } else {}
        if (class(list1) == "Enrichment" & class(list2) == "Enrichment") {
            temp1 <- list1
            temp2 <- list2
            if (nrow(temp1["Data"]) == 0) {
                warning('[Enrichment:compareEnrichment] "Enrichment" data is empty for list1', call. = FALSE)
                return(invisible())
            } else {}
            if (nrow(temp2["Data"]) == 0) {
                warning('[Enrichment:compareEnrich] "Enrichment" data is empty for list2', call. = FALSE)
                return(invisible())
            } else {}
        } else {
            if (class(list1) == "character" & class(list2) == "character") {
                cat("################ Init files ################\n")
                initFiles(pattern = pattern, snpInfoDir = snpInfoDir, signalFile = signalFile, ldThresh = ldThresh, LD = FALSE, mc.cores = mc.cores)
                cat("############### Read list 1 ################\n")
                temp1 <- readEnrichment(pattern = pattern, signalFile = signalFile, transcriptFile = transcriptFile, 
                                        snpListDir = list1, snpInfoDir = snpInfoDir, distThresh = distThresh, 
                                        sigThresh = sigThresh, LD = LD, extendMethod = extendMethod, mc.cores = mc.cores)
                reSample(object = temp1, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores)
                temp1 <- reset(temp1, "Resampling")
                cat("############### Read list 2 ################\n")
                temp2 <- readEnrichment(pattern = pattern, signalFile = signalFile, transcriptFile = transcriptFile, 
                                        snpListDir = list2, snpInfoDir = snpInfoDir, distThresh = distThresh, 
                                        sigThresh = sigThresh, LD = LD, extendMethod = extendMethod, mc.cores = mc.cores)
                reSample(object = temp2, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores)
                temp2 <- reset(temp2, "Resampling")
            } else {
                stop('[Enrichment:compareEnrichment] "Enrichment" object is required', call. = FALSE)
                return(invisible())
            }
        }
        l1 <- temp1@eSNP@List
        l2 <- temp2@eSNP@List
        if (identical(l1, l2)) {
            stop('[Enrichment:compareEnrichment] both lists are identical', call. = FALSE)
            return(invisible())
        } else {}
        if (length(l1)<length(l2)) {
            object1 <- temp1
            object2 <- temp2
            rm(temp1, temp2)
            namesRes <- c("Enrichment_1", "Enrichment_2", "PVALUE_1", "PVALUE_2", "PVALUE", "nSAMPLE", "SNP_1", "SNP_2")
        } else {
            object1 <- temp2
            object2 <- temp1
            rm(temp1, temp2)
            namesRes <- c("Enrichment_2", "Enrichment_1", "PVALUE_2", "PVALUE_1", "PVALUE", "nSAMPLE", "SNP_2", "SNP_1")
        }
        object2 <- reset(object2, "Data")
        object2 <- reset(object2, "LD")

        cat("############ Comparison Start ##############\n")
        listRes <- eval(parse(text = paste("list(", paste(paste("Chrom", seq(22), " = NULL",  sep = ""), collapse = ", "), ")", sep = "")))
        for (iChr in seq(22)) {
            cat("  **** Chromosome", if (nchar(iChr) == 1) {paste("0", iChr, sep = "")} else {iChr}, "Comparison Start ****\n")
            listRes[[iChr]] <- .compareEnrich(object1 = object1@Chromosomes[[iChr]], object2 = object2@Chromosomes[[iChr]], nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
            if (identical(sort(object1@Chromosomes[[iChr]]@eSNP@List), sort(object2@Chromosomes[[iChr]]@eSNP@List))) {
                listRes[[iChr]] <- reset(listRes[[iChr]], "Z")
                listRes[[iChr]] <- reset(listRes[[iChr]], "Resampling")
                listRes[[iChr]]@eSNP@PValue <- as.numeric(NA)
                listRes[[iChr]]@xSNP@PValue <- as.numeric(NA)
            } else {}
            cat("  ***** Chromosome", if (nchar(iChr) == 1) {paste("0", iChr, sep = "")} else {iChr}, "Comparison End *****\n")
        }

        cat("  ******* Genome Comparison Start  *******\n")
        result <- .compareEnrich(object1 = object1, object2 = object2, nSample = nSample, sigThresh = sigThresh, MAFpool = MAFpool, extendMethod = extendMethod, mc.cores = mc.cores)
        result@Chromosomes <- listRes
        cat("  ******** Genome Comparison End  ********\n")
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
    } else {
        stop('[Enrichment:compareEnrichment] "Enrichment" object is required', call. = FALSE)
        return(invisible())
    }
})
