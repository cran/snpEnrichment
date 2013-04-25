######################################################################
############################# Class SNP ##############################
############################## Creation ##############################
######################################################################


setClass(
    Class = "SNP", 
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
        Resampling = matrix(0, ncol = 0, nrow = 0)
    )
)


setMethod(f = "newSNP", signature = "ANY", definition = function(List, Table, EnrichmentRatio, Z, PValue, Resampling) {
    if (missing(List)) {List <- character()} else {}
    if (missing(Table)) {Table <- matrix(0, ncol = 2, nrow = 2)} else {}
    if (missing(EnrichmentRatio)) {EnrichmentRatio <- numeric()} else {}
    if (missing(Z)) {Z <- numeric()} else {}
    if (missing(PValue)) {PValue <- numeric()} else {}
    if (missing(Resampling)) {Resampling <- matrix(0, ncol = 5, nrow = 0)} else {}
    return(new("SNP", List = List, Table = Table, EnrichmentRatio = EnrichmentRatio, Z = Z, PValue = PValue, Resampling = Resampling))
})


setMethod(f = "is.SNP", signature = "ANY", definition = function(object){
    if (length(object)>1) {
        return(sapply(object, is.SNP))
    } else {
        if (class(object) == "SNP") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
})


.SNP.show <- function(object){
    cat("   - List :", paste("(", length(object@List), ")", sep = ""), if (length(object@List) <= 5) {if (length(object@List) == 0) {"NA"} else {object@List}} else {paste(paste(object@List[seq(5)], collapse = " "), "...")})
    cat("\n   - Table :", paste("(", paste(dim(object@Table), collapse = "x"), ")", sep = ""))
        if (all(object@Table == 0)) {cat(" NA")} else {print(object@Table, quote = FALSE)}
    cat("\n   - EnrichmentRatio :", ifelse(length(object@EnrichmentRatio) == 0, NA, object@EnrichmentRatio))
    cat("\n   - Z :", ifelse(length(object@Z) == 0, NA, object@Z))
    cat("\n   - PValue :", ifelse(length(object@PValue) == 0, NA, object@PValue))
    cat("\n   - Resampling :", paste("(", paste(dim(object@Resampling), collapse = "x"), ")", sep = ""), ifelse(nrow(object@Resampling) == 0, 0, nrow(object@Resampling)))
    cat("\n")
    return(invisible())
}
setMethod(f = "show", signature = "SNP", definition = function(object){cat("     ~~~ Class:", class(object), "~~~\n"); .SNP.show(object)})


setMethod(f = "[", signature = "SNP", definition = function(x, i, j, drop){
    switch(EXPR = i, 
        "List" = {return(x@List)}, 
        "Table" = {return(x@Table)}, 
        "EnrichmentRatio" = {return(x@EnrichmentRatio)}, 
        "Z" = {return(x@Z)}, 
        "PValue" = {return(x@PValue)}, 
        "Resampling" = {return(x@Resampling)},
        stop('[SNP:get] ', i, ' is not a "SNP" slot', call. = FALSE)
    )
    return(invisible())
})


setMethod(f = "[<-", signature = "SNP", definition = function(x, i, j, value){
    switch(EXPR = i, 
        "List" = {x@List <- value}, 
        "Table" = {x@Table <- value}, 
        "EnrichmentRatio" = {x@EnrichmentRatio <- value}, 
        "Z" = {x@Z <- value}, 
        "PValue" = {x@PValue <- value}, 
        "Resampling" = {x@Resampling <- value}, 
        stop('[SNP:set] ', i, ' is not a "SNP" slot', call. = FALSE)
    )
    validObject(x)
    return(x)
})
