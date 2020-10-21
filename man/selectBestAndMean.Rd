\name{selectBestAndMean}
\alias{selectBestAndMean}
\title{Select best HZEI and mean}
\usage{
selectBestAndMean(sequenceVector, clusterName, increaseHZEI=TRUE)
}
\arguments{
\item{sequenceVector}{Character vector of nucleotide sequences}
\item{clusterName}{Name of cluster generated with package parallel}
\item{increaseHZEI}{Logical value if HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
}
\description{
From all sequences of a generation report highest HZEI integral and mean HZEI integral of all.
}
\value{
Numeric vector of length 2 stating the best HZEI integral and the mean HZEI integral of a nucleotide sequence vector \code{sequenceVector}. Depending on the \code{increaseHZEI} mode, the best HZEI integral value is either the highest (for \code{increaseHZEI}==TRUE) or lowest (for \code{increaseHZEI}==FALSE).
}
\examples{
## Setup cluster
library(parallel)
nCores <- 1
clust <- makeCluster(nCores)
clusterExport(clust, list('getOverlappingVectorsFromVector', 'hex',
'calculateHZEIint'), envir = environment())
selectBestAndMean(c('CGCGATACGCTAAGCGCTACCGATAGTGGA','TGGGATATTTTAAGCGCTGACGATAGTGGA'),
clust, increaseHZEI=TRUE)
}
