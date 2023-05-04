#' Function for returning the signatures
#'
#' This function returns the genes in the signatures slot of the class. Optionally filtered by signature name
#' @title returnSignatures
#' @param x an object of class signatureTester
#' @param signatures a vector of different gene signatures to be tested (default NULL)
#' @keywords expression testing signatures
#' @export returnSignatures
#' @returns This returns a list of the signatures requested
#' @examples
#' returnSignatures(x = object)
#' returnSignatures(x = object, sigs = c("sig1", "sig2"))

returnSignatures <- function(x, Signatures = NULL) {
    if ( is.vector(Signatures)) {
        sigs = x@signatures[Signatures]
    } else {
        sigs = x@signatures
    }
    sigs
}

#' Function for filtering genes out of the signatures
#'
#' This function returns the signatures requested (default: all) and will also filter out the genes that are not in the stored eSet
#'
#' @title filterSignatureGenes
#' @param x an object of class signatureTester
#' @param Signatures a vector of different gene signatures to be tested (default NULL)
#' @param missingGenes what to do with the missing genes (currently "remove" or nothing)
#' @param verbose will output some metrics about the genes
#' @returns This returns a list of the signatures will the genes not found in the stored eSet removed.
#' @keywords expression testing signatures
#' @export filterSignatureGenes
#' @examples
#' filterSignatureGenes(x = object)
#' filterSignatureGenes(x = object, missingGenes = "remove", sigs = c("sig1", "sig2"))


filterSignatureGenes <- function(x, missingGenes = "remove", Signatures = NULL, verbose = FALSE) {
    ## get the signatures
    sigs <- returnSignatures(x, Signatures)
    ## iterate over the signatures and deal with genes missing from the expressionSet
    mapply(function(sig, name) {
        if (verbose) {
            original <- length(sig)
            remain <- sum(sig %in% featureNames(x@Expression))
            message("Signature: ", name, " - ", remain, " out of ", original, " genes in original dataset")
        }
        if (missingGenes == "remove") {
            sig[sig %in% featureNames(x@Expression)]
        } else {
            sig
        }
    }, sigs, names(sigs), SIMPLIFY = FALSE)
}

#' Function for removing missing genes from the signatures permanently
#'
#' This function returns a new copy of the class with the signatures filtered to only the genes that exist in the eSet
#'
#' This can be useful if you only want to store the genes that exist, but by default the downstream functions will soft filter these genes too.
#'
#' @title hardFilterSignatures
#' @param x an object of class signatureTester
#' @param Signatures a vector of different gene signatures to be tested (default NULL)
#' @param missingGenes what to do with the missing genes (currently "remove" or nothing)
#' @param verbose will output some metrics about the genes
#' @returns This returns a list of the signatures will the genes not found in the stored expressionSet removed.
#' @keywords expression testing signatures
#' @export hardFilterSignatures
#' @examples
#' hardFilterSignatures(x = object)
#' hardFilterSignatures(x = object, missingGenes = "remove", sigs = c("sig1", "sig2"))


hardFilterSignatures <- function(x, missingGenes = "remove", Signatures = NULL, verbose = FALSE) {
    x@signatures <- filterSignatureGenes(x,
                                         missingGenes = missingGenes,
                                         Signatures = Signatures,
                                         verbose = verbose)
    x@attributes$genefilter <- "Missing Genes Removed"
    x
}


















