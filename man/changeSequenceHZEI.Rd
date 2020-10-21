\name{changeSequenceHZEI}
\alias{changeSequenceHZEI}
\title{Adjust HZEI integral of nucleotide sequence}
\usage{
changeSequenceHZEI(inSeq, increaseHZEI=TRUE, nGenerations=30, parentSize=300,
startParentSize=1000, bestRate=50, semiLuckyRate=20, luckyRate=5, mutationRate=1e-04,
optiRate=100, sdMaximalHBS=10, acMaximalMaxent=4, nCores=-1)
}
\arguments{
\item{inSeq}{Character value of nucleotide sequence (min 24nt long, only bases A, G, T or C)}
\item{increaseHZEI}{Logical value if HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
\item{nGenerations}{Numeric value setting maximal number of generations}
\item{parentSize}{Numeric value setting size of parent generations, generated from previous generations}
\item{startParentSize}{Numeric value setting size of initiated parent generation of sequences}
\item{bestRate}{Numeric value setting percentage how many of the fittest sequences are used to produce the next generation}
\item{semiLuckyRate}{Numeric value setting percentage of sequences which are selected for breeding with a probability based on the respective HZEI-score integral}
\item{luckyRate}{Numeric value setting percentage of sequences which are randomly selected for breeding}
\item{mutationRate}{Numeric value setting chance of each codon, to mutate randomly within a child sequence}
\item{optiRate}{Numeric value setting level of HZEI integral optimization}
\item{sdMaximalHBS}{Numeric value of minimal HBS of SDs which should be tried to be degraded in their intrinsic strength}
\item{acMaximalMaxent}{Numeric value of minimal MaxEntScan score of SAs which should be tried to be degraded in their intrinsic strength}
\item{nCores}{Numeric value setting number of cores which should be used for parallel computations. If set to '-1' all availible cores are selected.}
}
\description{
Adjust the HZEI integral of a nucleotide sequence (min. 24nt long)
}
\value{
Character value of a nucleotide sequence encoding the same amino acid sequence as \code{inSeq}, but an increased HZEI integral, due to alternative codon selection.
}
\examples{
## Load R packages
library('parallel')
library('utils')
library('data.table')

## Set parameters for genetic algorithm
inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACC'

## Increase HZEI integral
res <- changeSequenceHZEI(inSeq)

## Setting additional parameters
res <- changeSequenceHZEI(inSeq, increaseHZEI=TRUE, nGenerations=30, parentSize=300,
startParentSize=1000, bestRate=50, semiLuckyRate=20, luckyRate=5, mutationRate=1e-04,
optiRate=100, sdMaximalHBS=10, acMaximalMaxent=4, nCores=1)

## Access sequence with highest generated HZEI intregral
res[[3]]

}
