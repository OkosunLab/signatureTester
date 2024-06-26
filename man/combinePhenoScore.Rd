% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreAnalysis.R
\name{combinePhenoScore}
\alias{combinePhenoScore}
\title{combinePhenoScore}
\usage{
combinePhenoScore(x, ...)
}
\arguments{
\item{x}{an object of class signatureTester}

\item{signatures}{a vector of the signatures to return (default = NULL)}
}
\value{
a dataframe of the phenodata from the eSet and the signature scores
}
\description{
Function for returning pheno data from the eSet along with the calculated scores
}
\details{
This function compares the scores and the groups
}
\examples{
combineExprsSurv(x = object)
}
\keyword{expression}
\keyword{signatures}
\keyword{testing}
