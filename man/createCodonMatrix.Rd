\name{createCodonMatrix}
\alias{createCodonMatrix}
\title{Create codon matrix from coding nucleotide sequence}
\usage{
createCodonMatrix(cds)
}
\arguments{
\item{cds}{Character value of nucleotide sequence whose HZEI integral will be calculated. It should be at least 3 nt long and only contain bases 'A', 'G', 'C', 'T'. Length must be a multiple of 3.}
}
\description{
This function creates a codon matrix with 2 rows and as many columns as codons within the sequence.
}
\value{
Character matrix holding the encoded codon sequence in both rows.
}
\examples{
## Example to create codon matrix
createCodonMatrix("ATGAATGATCAAAAGCTAGCC")

}
