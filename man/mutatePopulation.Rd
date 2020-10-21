\name{mutatePopulation}
\alias{mutatePopulation}
\title{Randomly exchange codons within a set of sequences.}
\usage{
mutatePopulation(sequenceVector, codonReplacementChance)
}
\arguments{
\item{sequenceVector}{Character vector of nucleotide sequences (at least 3 nt long)}
\item{codonReplacementChance}{Numeric value of chance of a codons within the sequences to get exchanged to another codon encoding the same amino acid}
}
\description{
For every codon within a set of nucleotide sequences randomly exchange the codon encoding the same amino acid to a certain chance.
}
\value{
Creates a character vector of coding nucleotide sequences encoding the same amino acid sequence as the entered \code{sequenceVector}. By a mutation rate stated in \code{codonReplacementChance}, codons are randomly exchanged, by alternative codons encoding the same amino acid.
}
\examples{
mutatePopulation(c("CGCGATACGCTAAGCGCTACCGATAGTGGA","TGGGATATTTTAAGCGCTGACGATAGTGGA"), 0.1)

}
