% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreAnalysis.R
\name{combineExprsSurv}
\alias{combineExprsSurv}
\title{combineExprsSurv}
\usage{
combineExprsSurv(x, Signatures = NULL, TIME = "OS_Time", IND = "OS_IND")
}
\arguments{
\item{x}{an object of class signatureTester}

\item{TIME}{the name of the time column in the phenodata (default = "OS_Time")}

\item{IND}{the name of the status indicator column in phenodata (default = "OS_IND")}

\item{signatures}{a vector of the signatures to return (default = NULL)}
}
\value{
The gene expression for the genes in the signatures along with the TIME and IND for survival analysis
}
\description{
Function for returning the gene expression along with survival stats
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
