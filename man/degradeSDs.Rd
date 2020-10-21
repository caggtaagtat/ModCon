\name{degradeSDs}
\alias{degradeSDs}
\title{Remove or degrade intrinisc strength of splice donors while adjusting HZEI integral.}
\usage{
degradeSDs(fanFunc, maxhbs=10, increaseHZEI=TRUE)
}
\arguments{
\item{fanFunc}{Codon matrix with two rows (see example below)}
\item{maxhbs}{Numeric treshold which strength of internal donor sites should be degraded}
\item{increaseHZEI}{Logical value of HZEI integral should be increased or decreased during SD degradation. If TRUE, function aims to increase HZEI integral.}
}
\description{
Degrade or remove splice donor sites of certain intrinsic strength (in HBS) from a coding sequence by codon selection.
}
\value{
Character value of a nucleotide sequence encoding the same amino acid as the entered  codon matrix \code{fanFunc}, but the intrinsic strength of all present splice donors (SD) sites is degraded as much as possible, in case they exceed the given treshold \code{maxhbs}.
}
\examples{
library(data.table)
## Initiaing the Codons matrix plus corresponding amino acids
ntSequence <- 'TTTTCGATCGGGATTAGCCTCCAGGTAAGTATCTATCGATCTATGCGATAG'
## Create Codon Matrix by splitting up the sequence by 3nt
fanFunc <- createCodonMatrix(ntSequence)
degradeSDs(fanFunc, maxhbs=10, increaseHZEI=TRUE)

}
