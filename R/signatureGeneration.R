#' Function for calculating a "mean expression" score
#'
#' This function will calculate the mean expression of the genes in each signature and store the resultant dataframe in the scores slot of the class
#' @title calculateMeanExpression
#' @param x an object of class signatureTester
#' @param signatures a vector of different gene signatures to be tested (default NULL)
#' @keywords expression testing signatures
#' @export calculateMeanExpression
#' @returns This a copy of the class with the mean expression of each signature in the scores slot
#' @examples
#' calculateMeanExpression(x = object)
#' calculateMeanExpression(x = object, sigs = c("sig1", "sig2"))

calculateMeanExpression <- function(x, Signatures = NULL) {
    ## get the signatures
    if ( is.null(x@attributes$genefilter) ) {
        warning("Some genes may be in the signatures slot that are not in the cohort.\n",
                "These will be ignored when generating the signature.\n")
        sigs <- filterSignatureGenes(x, Signatures = Signatures, verbose = TRUE)
    } else {
        sigs <- returnSignatures(x, Signatures)
    }
    ## calculate mean expression of the genes
    x@scores <- sapply(sigs, function(genes) {
        colMeans(exprs(x@Expression[genes,]))
    }) %>% as.data.frame()
    x
}
