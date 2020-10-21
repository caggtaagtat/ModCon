\name{generateRandomCodonsPerAA}
\alias{generateRandomCodonsPerAA}
\title{Randomly choose Codon to encode amino acid sequence}
\usage{
generateRandomCodonsPerAA(aaVector)
}
\arguments{
\item{aaVector}{Character vector of amino acids in three lettre code (e.g. Met)}
}
\description{
Encode amino acid sequence by random codon selection
}
\value{
Character value of a nucleotide sequence encoding the same amino acid as the entered by \code{aaVector} by random Codon selection.
}
\examples{
generateRandomCodonsPerAA(c('Lys','Lys'))

}
