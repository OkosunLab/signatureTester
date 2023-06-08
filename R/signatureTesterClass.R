#' Class for testing expression signatures
#'
#' This class will hold your segments and sample information
#' @title signatureTester
#' @slot rawExpression this will store the initial expression set that is built when the class is created.
#' @slot signatures this is a list of the signatures to be tested.
#' @keywords expression testing signatures
#' @export signatureTester
#' @examples
#' signatureTester()

signatureTester <- setClass("signatureTester",
                            slots = c(
                                Expression="ExpressionSet",
                                scores="data.frame",
                                signatures="list",
                                attributes = "list",
                                assignments = "list"
                            ))


setMethod("show", "signatureTester",
          function(object) {

              cat("Expression Data:\n\tfeatures:",
                  nrow(object@Expression),
                  "\n\tsamples:",
                  ncol(object@Expression),
                  "\nSignatures:",
                  ifelse(!is.null(object@signatures),
                         cat(
                             paste("\n\tSignatures:", length(object@signatures)),
                             paste("\n\t\t", names(object@signatures))),
                         "\n\tno signatures stored")
                  ifelse(!is.null(object@scores),
                         cat(
                             paste("\n\tScores:", length(object@scores)),
                             paste("\n\t\t", names(object@scores))),
                         "\n\tno signatures stored"))
              }
              )

#' Function for building the class from an eSet
#'
#' This function will generate an instance of class "signatureTester"
#' @title signatureTesterFromExpressionSet
#' @param eSet this is an expressionset containing the gene expression of the data you are testing the signature in
#' @param signatures a list of different gene signatures to be tested (default NULL)
#' @keywords expression testing signatures
#' @export signatureTesterFromExpressionSet
#' @examples
#' signatureTesterFromExpressionSet(eSet = eSet)
#' signatureTesterFromExpressionSet(eSet = eSet, sigs = Signatures)


signatureTesterFromExpressionSet <- function(eSet, sigs = NULL) {
    if (is.null(sigs)) {
        signatureTester(Expression = eSet)
    } else if (is.data.frame(sigs)) {
        message("Signatures are dataframe, these will be converted to a list")
        sigs <- sigs %>% as.list()
        signatureTester(Expression = eSet,
                        signatures = sigs)
    } else if (is.list(sigs)) {
        signatureTester(Expression = eSet,
                        signatures = sigs)
    }
}
