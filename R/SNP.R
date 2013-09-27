######################################################################
########################### Class EnrichSNP ##########################
############################## Creation ##############################
######################################################################


setClass(
    Class = "EnrichSNP", 
    representation = representation(
        List = "character", 
        Table = "matrix", 
        EnrichmentRatio = "numeric", 
        Z = "numeric", 
        PValue = "numeric", 
        Resampling = "matrix"
    ), 
    prototype = prototype(
        List = character(), 
        Table = matrix(0, ncol = 2, nrow = 2), 
        EnrichmentRatio = numeric(), 
        Z = numeric(), 
        PValue = numeric(), 
        Resampling = matrix(0, ncol = 5, nrow = 0)
    )
)


setMethod(f = "enrichSNP", signature = "ANY", definition = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) {
    if (missing(List)) {
        List <- character()
    } else {}
    if (missing(Table)) {
        Table <- matrix(0, ncol = 2, nrow = 2)
    } else {}
    if (missing(EnrichmentRatio)) {
        EnrichmentRatio <- numeric()
    } else {}
    if (missing(Z)) {
        Z <- numeric()
    } else {}
    if (missing(PValue)) {
        PValue <- numeric()
    } else {}
    if (missing(Resampling)) {
        Resampling <- matrix(0, ncol = 5, nrow = 0)
    } else {}
    return(new("EnrichSNP", List = List, Table = Table, EnrichmentRatio = EnrichmentRatio, Z = Z, PValue = PValue, Resampling = Resampling))
})


setMethod(f = "is.EnrichSNP", signature = "ANY", definition = function(object){
    if (length(object)>1) {
        return(sapply(object, is.EnrichSNP))
    } else {
        if (class(object) == "EnrichSNP") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
})


.EnrichSNP.show <- function(object){
    cat("   - List :", 
        paste0("(", length(object@List), ")"), 
        if (length(object@List) <= 5) {
            if (length(object@List) == 0) {
                "NA"
            } else {
                object@List
            }
        } else {
            paste(paste(object@List[seq(5)], collapse = " "), "...")
        }
    )
    cat("\n   - Table :", paste0("(", paste(dim(object@Table), collapse = "x"), ")"))
        if (all(object@Table == 0)) {
            cat(" NA")
        } else {
            print(object@Table, quote = FALSE)
        }
    cat("\n   - EnrichmentRatio :", ifelse(length(object@EnrichmentRatio) == 0, NA, object@EnrichmentRatio))
    cat("\n   - Z               :", ifelse(length(object@Z) == 0, NA, object@Z))
    cat("\n   - PValue          :", ifelse(length(object@PValue) == 0, NA, object@PValue))
    cat("\n   - Resampling      :", paste0("(", paste(dim(object@Resampling), collapse = "x"), ")"), ifelse(nrow(object@Resampling) == 0, 0, nrow(object@Resampling)))
    cat("\n")
    return(invisible())
}
setMethod(f = "show", signature = "EnrichSNP", definition = function(object){cat("     ~~~ Class:", class(object), "~~~\n"); .EnrichSNP.show(object)})
setMethod(f = "print", signature = "EnrichSNP", definition = function(x, ...){cat("     ~~~ Class:", class(x), "~~~\n"); .EnrichSNP.show(x)})


setMethod(f = "[", signature = "EnrichSNP", definition = function(x, i, j, drop){
    switch(EXPR = i, 
        "List" = {return(x@List)}, 
        "Table" = {return(x@Table)}, 
        "EnrichmentRatio" = {return(x@EnrichmentRatio)}, 
        "Z" = {return(x@Z)}, 
        "PValue" = {return(x@PValue)}, 
        "Resampling" = {return(x@Resampling)}, 
        stop('[EnrichSNP:get] ', i, ' is not a "EnrichSNP" slot', call. = FALSE)
    )
    return(invisible())
})


setMethod(f = "[<-", signature = "EnrichSNP", definition = function(x, i, j, value){
    switch(EXPR = i, 
        "List" = {x@List <- value}, 
        "Table" = {x@Table <- value}, 
        "EnrichmentRatio" = {x@EnrichmentRatio <- value}, 
        "Z" = {x@Z <- value}, 
        "PValue" = {x@PValue <- value}, 
        "Resampling" = {x@Resampling <- value}, 
        stop('[EnrichSNP:set] ', i, ' is not a "EnrichSNP" slot', call. = FALSE)
    )
    validObject(x)
    return(x)
})

setMethod(f = "reset", signature = "EnrichSNP", definition = function(object, i) {
    switch(EXPR = i, 
        "List" = {object@List <- character()}, 
        "Table" = {object@Table <- matrix(0, ncol = 2, nrow = 2)}, 
        "EnrichmentRatio" = {object@EnrichmentRatio <- numeric()}, 
        "Z" = {object@Z <- numeric()}, 
        "PValue" = {object@PValue <- numeric()}, 
        "Resampling" = {object@Resampling <- matrix(0, ncol = 5, nrow = 0)}, 
        stop('[EnrichSNP:reset] ', i, ' is not a "EnrichSNP" slot', call. = FALSE)
    )
    return(object)
})
