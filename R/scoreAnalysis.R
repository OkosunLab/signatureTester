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

thresholdScores <- function(x, method = "ntile", n = 2) {
    if (method == "ntile") {
        x@assignments$groups <- apply(x@scores, 2, function(col) {
            ntile(col, n)}) %>%
            as.data.frame() %>%
            mutate(Mixture = rownames(x@scores)) %>%
            column_to_rownames("Mixture")
    } else if (method == "mean") {

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
