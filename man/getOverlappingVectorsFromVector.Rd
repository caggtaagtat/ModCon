\name{getOverlappingVectorsFromVector}
\alias{getOverlappingVectorsFromVector}
\title{Create overlapping subvectors}
\usage{
getOverlappingVectorsFromVector(largeVector, subvectorLength, subvectorOverlap )
}
\arguments{
\item{largeVector}{Large character vector to break down into overlapping subvectors}
\item{subvectorLength}{Numeric value of length of smaller subvectors}
\item{subvectorOverlap}{Numeric value of length of subvector overlap}
}
\description{
Create overlapping subvectors from large vector
}
\value{
Creates a list of overlapping subvectors from an input vector \code{largeVector}. The length of these overlapping subvectors is stated by \code{subvectorLength} and the overlap of the resulting subvectors is stated by \code{subvectorOverlap}.
}
\examples{
getOverlappingVectorsFromVector(c(1,2,3,4), 2, 1)

}
