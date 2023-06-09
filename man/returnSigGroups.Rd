% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreAnalysis.R
\name{returnSigGroups}
\alias{returnSigGroups}
\title{returnSigGroups}
\usage{
returnSigGroups(x, Signatures = NULL)
}
\arguments{
\item{x}{an object of class signatureTester}

\item{signatures}{a vector of different gene signatures to be tested (default NULL)}
}
\value{
a dataframe of the assigned groups
}
\description{
Function for returning the assigned groups for signatures
}
\details{
This function compares the scores and the groups
}
\examples{
returnSigGroups(x = object)
}
\keyword{expression}
\keyword{signatures}
\keyword{testing}
