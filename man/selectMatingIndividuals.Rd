\name{selectMatingIndividuals}
\alias{selectMatingIndividuals}
\title{Selecting mating sequences from a pool of sequences}
\usage{
selectMatingIndividuals(inputGeneration, whoMatesBestPercent=40, whoMatesSemiRandom=20,
whoMatesLuckily=5, clust, increaseHZEI=TRUE)
}
\arguments{
\item{inputGeneration}{Character vector of nucleotide sequences}
\item{whoMatesBestPercent}{Numeric value e.g. 20 (which would mean that sequences with the top 20 percent highest HZEI integral are selected for mating)}
\item{whoMatesSemiRandom}{Numeric value (is always lower than total number of sequences in input_generation)}
\item{whoMatesLuckily}{Numeric value (is always lower than total number of sequences in input_generation)}
\item{clust}{Name of cluster generated with package parallel}
\item{increaseHZEI}{Logical value of HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
}
\description{
Selecting sequences from a pool of nucleotide sequences based in chance and their HZEI integral.
}
\value{
Character vector of nucleotide sequences which are selected from an entered vector of nucleotide sequences \code{inputGeneration} for creation of filial sequences by recombination. Sequences are selected by different criteria stated by  \code{whoMatesBestPercent}, \code{whoMatesSemiRandom}, \code{whoMatesLuckily} and \code{increaseHZEI}.
}
\examples{
## Setup cluster
library(parallel)
nCores <- 1
clust <- makeCluster(nCores)
clusterExport(clust, list('getOverlappingVectorsFromVector',
'hex'), envir=environment())
selectMatingIndividuals(c('CGCGATACGCGCGATACG','CGCGATACGTGGGATATT',
'CTAAGCGCTCGCGATACG','CGCGATACGTTAAGCGCT','GACGATAGTCGCGATACG'),
 whoMatesBestPercent=40, whoMatesSemiRandom=1, whoMatesLuckily=1, clust, increaseHZEI=TRUE)

}
