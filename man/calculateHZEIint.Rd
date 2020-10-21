\name{calculateHZEIint}
\alias{calculateHZEIint}
\title{Calculate HZEI integral of nucleotide sequence}
\usage{
calculateHZEIint(ntSequence)
}
\arguments{
\item{ntSequence}{Character value of nucleotide sequence whose HZEI integral
will be calculated. It should be at least 11 nt long and only contain bases
'A', 'G', 'C', 'T'.}
}
\description{
This function calcuales the HZEI integral of a nucleotide sequence.
}
\value{
Integer value stating the HZEI integral of the given sequence \code{ntSequence}
}
\examples{
## Example to increase HZEI integral for a given coding sequence
x <- calculateHZEIint('ATACCAGCCAGCTATTACATTT')

}
