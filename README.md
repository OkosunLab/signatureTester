# signatureTester

This package has been built to simply the testing of gene signatures in bulk expression data.

## Quickstart

To install the package make sure devtools in installed

```r
if (! require(devtools)) {
  install.packages("devtools")
  }
```

Then install the package directly from github:

```r
devtools::install_github("OkosunLab/signatureTester")
```

### Testing your first signatures

You need two things to get started:
  - Expression data stored as an [expression set](bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/BiobaseDevelopment.html)
  - Your gene signatures stored as a named list or a dataframe

A simplified workflow might look like this:

```r
## load eSet
expressionSet <- readRDS("eSet.rds")
## Generate the signatures 
signatures <- list(
  "sig1" = c("GeneA", "GeneB"),
  "sig2" = c("GeneC", "GeneD")
  )
## Build the class 
SignatureTester <- signatureTesterFromExpressionSet(eSet = expressionSet, sigs = signatures)
## Calculate the mean expression of all the genes in these signatures
SignatureTester <- calculateMeanExpression(SignatureTester)
## Threshold the scores into two groups, using the ntile method (i.e. split at the median value)
thresholdScores(method = "ntile", n = 2)
```
