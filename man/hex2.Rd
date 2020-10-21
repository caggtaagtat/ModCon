\docType{data}
\name{hex2}
\alias{hex2}
\title{Hexamers and Z scores for every 9-mer nucleotide sequence}
\format{A data frame with columns:
\describe{
 \item{seq}{Sequence of first hexamer (first 2 codons).}
 \item{seq2}{Sequence of last hexamer  (last 2 codons).}
 \item{first}{First codon within the hexamer.}
 \item{first_AA}{First encoded amino acid within the hexamer (three lettre code).}
 \item{second}{Second codon within the hexamer.}
 \item{second_AA}{Second encoded amino acid within the hexamer (three lettre code).}
 \item{third}{Second codon within the hexamer.}
 \item{third_AA}{Second encoded amino acid within the hexamer (three lettre code).}
 \item{hex1}{ZEI-score of first hexamer from HEXplorer.}
 \item{hex2}{ZEI-score of last hexamer from HEXplorer.}
 \item{sum_hex}{Sum of ZEI-scores of every overlapping hexamer}
 \item{code_AA}{All 3 encoded amino acid within the 9-mer nucleotid sequence.}
 \item{ss}{9-mer nucleotide sequence.}
 \item{gc}{Number of G or C nucleotides in the 9-mer sequence.}
}}
\usage{
hex2
}
\description{
Hexamers and Z scores for every 9-mer nucleotide sequence
}
\examples{
 hex2
}
\keyword{datasets}
