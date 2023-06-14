#' Function for zscoring scores
#'
#' This function will zscore the individual signature scores and return a new copy of the class.
#' @title scaleScores
#' @param x an object of class signatureTester
#' @keywords expression testing signatures
#' @export scaleScores
#' @returns This returns the class with the scores slot z scored
#' @examples
#' scaleScores(x = object)

scaleScores <- function(x) {
    x@scores <- scale(x@scores) %>% as.data.frame()
    x@attributes$scale = "zscored"
    x
}

#' Function for thresholding scores
#'
#' This function offers a few ways of thresholding the scores found in the scores slot of the class
#' @title thresholdScores
#' @param x an object of class signatureTester
#' @param method the method for thresholding (current options: "ntile")
#' @param n the number of groups, if relevant (default: 2)
#' @keywords expression testing signatures
#' @export thresholdScores
#' @returns This returns the class with groups added to the assignment slot
#' @examples
#' thresholdScores(x = object, method = "ntile", n = 2)

thresholdScores <- function(x, method = "ntile", n = 2, thres = NULL) {
    if (method == "ntile") {
        x@assignments$groups <- apply(x@scores, 2, function(col) {
            ntile(col, n)}) %>%
            as.data.frame() %>%
            mutate(Mixture = rownames(x@scores)) %>%
            column_to_rownames("Mixture")
    } else if (method == "mean") {
        x@assignments$groups <- apply(x@scores, 2, function(col) {
            ifelse(col > mean(col), "high" , "low" ) %>%
            as.data.frame() %>%
            mutate(Mixture = rownames(x@scores)) %>%
            column_to_rownames("Mixture")
    })
    }
    x
}

#' Function for combining the scores and the assigned groups
#'
#' This function grabs the scores and the groups and creates a dataset so they can be compared.
#' @title combineGroupsScores
#' @param x an object of class signatureTester
#' @keywords expression testing signatures
#' @export combineGroupsScores
#' @returns a dataframe of the scores and groups combined.
#' @examples
#' combineGroupsScores(x = object)

combineGroupsScores <- function(x) {
    left_join(x@assignments$groups %>%
                  rownames_to_column("ID") %>%
                  pivot_longer(-ID, names_to = "Signatures", values_to = "Group"),
              x@scores %>%
                  rownames_to_column("ID") %>%
                  pivot_longer(-ID, names_to = "Signatures", values_to = "Score"))
}

#' Function for comparing the scores and the assigned groups
#'
#' This function compares the scores and the groups
#' @title compareGroupScores
#' @param x an object of class signatureTester
#' @keywords expression testing signatures
#' @export compareGroupScores
#' @returns a boxplot of the scores across the groups in each signature
#' @examples
#' compareGroupScores(x = object)

compareGroupScores <- function(x) {
    combineGroupsScores(x) %>%
        ggplot(aes(x = Group, y = Score, colour = as.character(Group))) +
        geom_boxplot() +
        facet_wrap(~Signatures, scales = "free_y")
}


#' Function for returning the assigned groups for signatures
#'
#' This function compares the scores and the groups
#' @title returnSigGroups
#' @param x an object of class signatureTester
#' @param signatures a vector of different gene signatures to be tested (default NULL)
#' @keywords expression testing signatures
#' @export compareGroupScores
#' @returns a dataframe of the assigned groups
#' @examples
#' returnSigGroups(x = object)

returnSigGroups <- function(x, Signatures = NULL) {
    if ( is.null(Signatures) ) {
        x@assignments$groups
    } else {
        select(x@assignments$groups, Signatures)
    }
}

#' Function for returning columns from the expression set phenotype data
#'
#' This function compares the scores and the groups
#' @title compareGroupScores
#' @param x an object of class signatureTester
#' @param columns a vector of the columns to return (default = NULL)
#' @keywords expression testing signatures
#' @export compareGroupScores
#' @returns a data frame of the specified columns from the pData of the eSet
#' @examples
#' returnPhenoColumns(x = object)

returnPhenoColumns <- function(x, Columns = NULL) {
    if ( is.null(Columns) ) {
        pData(x@Expression)
    } else {
        pData(x@Expression)[,Columns]
    }
}

#' Function for returning the gene expression along with assigned groups
#'
#' This function compares the scores and the groups
#' @title combineExprsGroup
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @keywords expression testing signatures
#' @export combineExprsGroup
#' @returns a data frame of the expression values for the signature along with the identified groups
#' @examples
#' combineExprsGroup(x = object)

combineExprsGroup <- function(x, Signatures = NULL) {
    left_join(t(returnSigGeneExpression(x, Signatures = Signatures)) %>% as.data.frame() %>% rownames_to_column("ID"),
              returnSigGroups(x, Signatures = Signatures) %>%  rownames_to_column("ID"))
}

#' Function for returning the genes in the specified signatures
#'
#' This function compares the scores and the groups
#' @title returnSignatureGenes
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @keywords expression testing signatures
#' @export returnSignatureGenes
#' @returns a list of the genes in the specified signatures
#' @examples
#' returnSignatureGenes(x = object)

returnSignatureGenes <- function(x, Signatures = NULL) {
    filterSignatureGenes(x, Signatures = Signatures) %>%
        unlist() %>%
        unique()
}

#' Function for returning the gene expression of genes in the specified signatures
#'
#' This function compares the scores and the groups
#' @title returnSigGeneExpression
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @keywords expression testing signatures
#' @export returnSigGeneExpression
#' @returns The gene expression values for the genes in the specified signatures
#' @examples
#' returnSigGeneExpression(x = object)

returnSigGeneExpression <- function(x, Signatures = NULL) {
    Genes <- returnSignatureGenes(x, Signatures = Signatures)
    exprs(x@Expression[Genes,])
}

#' Function for returning the gene expression along with survival stats
#'
#' This function compares the scores and the groups
#' @title combineExprsSurv
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @param TIME the name of the time column in the phenodata (default = "OS_Time")
#' @param IND the name of the status indicator column in phenodata (default = "OS_IND")
#' @keywords expression testing signatures
#' @export combineExprsSurv
#' @returns The gene expression for the genes in the signatures along with the TIME and IND for survival analysis
#' @examples
#' combineExprsSurv(x = object)

combineExprsSurv <- function(x, Signatures = NULL, TIME = "OS_Time", IND = "OS_IND") {
    left_join(t(returnSigGeneExpression(x, Signatures = Signatures)) %>% as.data.frame() %>% rownames_to_column("ID"),
              returnPhenoColumns(x, Columns = c(TIME, IND)) %>% rownames_to_column("ID"))
}

#' Function for returning pheno data from the eSet along with the assigned groups
#'
#' This function compares the scores and the groups
#' @title combinePhenoGroup
#' @param x an object of class signatureTester
#' @keywords expression testing signatures
#' @export combinePhenoGroup
#' @returns a dataframe of the phenodata from the eSet and the signature groups
#' @examples
#' combinePhenoGroup(x = object)

combinePhenoGroup <- function(x) {
    left_join(pData(x@Expression) %>% rownames_to_column("ID"),
              x@assignments$groups %>% rownames_to_column("ID"))
}

#' Function for returning pheno data from the eSet along with the calculated scores
#'
#' This function compares the scores and the groups
#' @title combinePhenoScore
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @keywords expression testing signatures
#' @export combinePhenoScore
#' @returns a dataframe of the phenodata from the eSet and the signature scores
#' @examples
#' combineExprsSurv(x = object)

combinePhenoScore <- function(x) {
    left_join(pData(x@Expression) %>% rownames_to_column("ID"),
              x@scores %>% rownames_to_column("ID"))
}

#' Function for returning the gene expression along with survival stats
#'
#' This function compares the scores and the groups
#' @title calculateRatio
#' @param x an object of class signatureTester
#' @param signatures a vector of the signatures to return (default = NULL)
#' @param sig1 the name of the signature to be the numerator
#' @param sig2 the name of the signature to be the denominator
#' @param name the name of the new score that will be added to the scores slot (by default this will be sig1.over.sig2)
#' @keywords expression testing signatures
#' @export calculateRatio
#' @returns a boxplot of the scores across the groups in each signature
#' @examples
#' calculateRatio(x = object, sig1 = "Signature1", sig2 = "Signature2")
#' calculateRatio(x = object, sig1 = "Signature1", sig2 = "Signature2", name = "SigRatio")

calculateRatio <- function(x, sig1, sig2, name = NULL) {
    ratio <- x@scores[sig1]/x@scores[sig2]
    if ( is.null(name) ) {
        x@scores[paste(sig1, "over", sig2, sep = ".")] <- ratio
    } else {
        x@scores[name] <- ratio
    }
    x
}
