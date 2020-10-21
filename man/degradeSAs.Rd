\name{degradeSAs}
\alias{degradeSAs}
\title{Remove or degrade intrinisc strength of splice acceptors while adjusting HZEI integral.}
\usage{
degradeSAs(fanFunc, maxhbs=10, maxME=4, increaseHZEI=TRUE)
}
\arguments{
\item{fanFunc}{codon matrix with two rows (see example below)}
\item{maxhbs}{Numeric treshold which strength of internal donor sites should be degraded (in HBS)}
\item{maxME}{Numeric treshold which strength of internal acceptor sites should be degraded (in MaxEntScan score)}
\item{increaseHZEI}{Logical value if HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
}
\description{
Degrade or remove splice acceptor sites of certain intrinsic strength (in MaxEntScan score) from a coding sequence by codon selection while keeping the HZEI integral up.
}
\value{
Character value of a nucleotide sequence encoding the same amino acid as the entered  codon matrix \code{fan}, but the intrinsic strength of all present splice acceptor (SA) sites is degraded as much as possible, in case they exceed the given treshold \code{maxME}. Additionally, splice donor site strengths greater \code{maxhbs} are avoided, during SA degradation.
}
\examples{
library(data.table)
sdMaximalHBS <- 10
acMaximalMaxent <- 4
increaseHZEI <- TRUE
## Initiaing the Codons matrix plus corresponding amino acids
ntSequence <- 'TTTTGTCTTTTTCTGTGTGGCAGTGGGATTAGCCTCCTATCGATCTATGCGATA'
## Create Codon Matrix by splitting up the sequence by 3nt
fanFunc <- createCodonMatrix(ntSequence)
degradeSAs(fanFunc, maxhbs=sdMaximalHBS, maxME=acMaximalMaxent, increaseHZEI=increaseHZEI)

}
