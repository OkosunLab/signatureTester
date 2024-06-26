% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreAnalysis.R
\name{calculateRatio}
\alias{calculateRatio}
\title{calculateRatio}
\usage{
calculateRatio(x, sig1, sig2, name = NULL)
}
\arguments{
\item{x}{an object of class signatureTester}

\item{sig1}{the name of the signature to be the numerator}

\item{sig2}{the name of the signature to be the denominator}

\item{name}{the name of the new score that will be added to the scores slot (by default this will be sig1.over.sig2)}

\item{signatures}{a vector of the signatures to return (default = NULL)}
}
\value{
a boxplot of the scores across the groups in each signature
}
\description{
Function for generating a ration of two scores
}
\details{
This function compares the scores and the groups
}
\examples{
calculateRatio(x = object, sig1 = "Signature1", sig2 = "Signature2")
calculateRatio(x = object, sig1 = "Signature1", sig2 = "Signature2", name = "SigRatio")
}
\keyword{expression}
\keyword{signatures}
\keyword{testing}
