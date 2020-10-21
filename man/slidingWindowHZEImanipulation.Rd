\name{slidingWindowHZEImanipulation}
\alias{slidingWindowHZEImanipulation}
\title{Quickly manipulate HZEI integral of nucleotide sequence}
\usage{
slidingWindowHZEImanipulation(inSeq, increaseHZEI=TRUE)
}
\arguments{
\item{inSeq}{Character value of nucleotide sequence (min 21nt long, only bases 'A', 'G', 'T' or 'C')}
\item{increaseHZEI}{Logical value if HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
}
\description{
Quickly manipulate HZEI integral of nucleotide sequence (min. 21nt long)
}
\value{
Character value of a nucleotide sequence encoding the same amino acid sequence as \code{inSeq}, but an increased HZEI integral, due to alternative codon selection, accomplished through sliding window optimization.
}
\examples{
# Load R packages
library('parallel')
library('utils')
library('data.table')

# Set parameters for genetic algorithm
inSeq <- 'ATGGAAGACGCCAAAAACATAAAGAAAGGCAGGCTAAGCCTAGCTTGCCATTGCCCGGCGCCATTCTATCCGCTGGAAGATGGAATT'

maximizedHZEIseq <- slidingWindowHZEImanipulation(inSeq, increaseHZEI=TRUE)
minimizedHZEIseq <- slidingWindowHZEImanipulation(inSeq, increaseHZEI=FALSE)

#Access sequence with maximized HZEI intregral
maximizedHZEIseq

#Access sequence with minimized HZEI intregral
minimizedHZEIseq

}
