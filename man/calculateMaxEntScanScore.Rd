\name{calculateMaxEntScanScore}
\alias{calculateMaxEntScanScore}
\title{Calculate MaxEntScan score of a splice site sequence}
\usage{
calculateMaxEntScanScore(seqVector, ssType)
}
\arguments{
\item{seqVector}{Character value of nucleotide sequence of a splice site sequence. SA sequences should be 23nt long (20 intronic, 3 exonic) and SD sequences should be 9nt long (3 exonic, 6 intronic). Only bases 'A', 'G', 'C', 'T' permitted.}
\item{ssType}{Numeric value which indicates the type of splice site. Either '3' for an SA or '5' for an SD.}
}
\description{
This function calculates the MaxEntScan score of either splice donor (SD) or acceptor sequences (SA).
}
\value{
Numeric vector stating the MaxEntScan score per splice site sequence entered with \code{seqVector}
}
\examples{
calculateMaxEntScanScore('TTCCAAACGAACTTTTGTAGGGA',3)
calculateMaxEntScanScore('GAGGTAAGT',5)

}
