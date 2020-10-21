\name{ModCon}
\alias{ModCon}
\title{ModCon}
\usage{
ModCon(cds, sdSeqStartPosition, upChangeCodonsIn=16, downChangeCodonsIn=16,
modconMode="optimalContext", sdMaximalHBS=10, acMaximalMaxent=4, optiRate=100,
nGenerations=30, parentSize=300, startParentSize=1000, bestRate=40,
semiLuckyRate=20, luckyRate=5, mutationRate=1e-04, nCores=-1)
}
\arguments{
\item{cds}{Character value of coding nucleotide sequence which holds the splice site of interest}
\item{sdSeqStartPosition}{Numeric value of the position of the first nucleotide of the splice donor of interest}
\item{upChangeCodonsIn}{Numeric value of number of codons to change upstream of the donor site of interest}
\item{downChangeCodonsIn}{Numeric value of number of codons to change downstream of the donor site of interest}
\item{modconMode}{Character value which determines, if 'optimalContext' (the default) the donor context will be adjusted to increase the HZEI weight, if changed to 'suboptimalContext', the HZEI weight will be decreased.}
\item{sdMaximalHBS}{Numeric value of minimal HBS of SDs which should be tried to be degraded in their intrinsic strength}
\item{acMaximalMaxent}{Numeric value of minimal MaxEntScan score of SAs which should be tried to be degraded in their intrinsic strength}
\item{optiRate}{Numeric value setting level of HZEI integral optimization}
\item{nGenerations}{Numeric value setting maximal number of generations}
\item{parentSize}{Numeric value setting size of parent generations, generated from previous generations}
\item{startParentSize}{Numeric value setting size of initiated parent generation of sequences}
\item{bestRate}{Numeric value setting percentage how many of the fittest sequences are used to produce the next generation}
\item{semiLuckyRate}{Numeric value setting percentage of sequences which are selected for breeding with a probability based on the respective HZEI-score integral}
\item{luckyRate}{Numeric value setting percentage of sequences which are randomly selected for breeding}
\item{mutationRate}{Numeric value setting chance of each codon, to mutate randomly within a child sequence}
\item{nCores}{Numeric value setting number of cores which should be used for parallel computations. If set to '-1' all availible cores are selected.}

}
\description{
Execute ModCon on a donor site within a coding sequnece either increasing or decreasing its HZEI weight.
}
\value{
Creates a character value of a coding nucleotide sequence encoding the same amino acid sequence as the entered \code{cds}, but with an alternative nucleotide surrounding around the splice donor (SD) sequence position, stated with \code{sdSeqStartPosition}. Depending on the entered \code{modconMode}, the SD surrounding is either adjusted aiming to enhance or decrease the splice site HEXplorer wheigth.
}
\examples{
## Load R packages
library('parallel')
library('utils')
library('data.table')

## Set parameters for simplest use of ModCon (optimizing to 100\%)
cds <- paste0('ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACCGCTGGAGAGCAACTGCA',

'TAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATTGCTTTTACAGATGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAA',
'TGTCCGTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGAATCGTCGTATGCAGTGAAAACTCTCTTCAATTCTTTAT',
'GCCGGTGTTGGGCGCGTTATTTATCGGAGTTGCAGTTGCGCCCGCGAACGACATTTATAATGAACGTGAATTGCTCAACAGTATGGGCATTTCG',
'CAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAAAAAAAGCTCCCAATCATCCAAAAAATTATTATCATGG',
'ATTCTAAAACGGATTACCAGGGATTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATACGATTTTGTGCCAGA',
'GTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGATCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACTGCC',
'TGCGTGAGATTCTCGCATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTTGTTCCATTCCATCACGGTT',
'TTGGAATGTTTACTACACTCGGATATTTGATATGTGGATTTCGAGTCGTCTTAATGTATAGATTTGAAGAAGAGCTGTTTCTGAGGAGCCTTCA',
'GGATTACAAGATTCAAAGTGCGCTGCTGGTGCCAACCCTATTCTCCTTCTTCGCCAAAAGCACTCTGATTGACAAATACGATTTATCTAATTTA',
'CACGAAATTGCTTCTGGTGGCGCTCCCCTCTCTAAGGAAGTCGGGGAAGCGGTTGCCAAGAGGTTCCATCTGCCAGGTATCAGGCAAGGATATG',
'GGCTCACTGAGACTACATCAGCTATTCTGATTACACCCGAGGGGGATGATAAACCGGGCGCGGTCGGTAAAGTTGTTCCATTTTTTGAAGCGAA',
'GGTTGTGGATCTGGATACCGGGAAAACGCTGGGCGTTAATCAAAGAGGCGAACTGTGTGTGAGAGGTCCTATGATTATGTCCGGTTATGTAAAC',
'AATCCGGAAGCGACCAACGCCTTGATTGACAAGGATGGATGGCTACATTCTGGAGACATAGCTTACTGGGACGAAGACGAACACTTCTTCATCG',
'TTGACCGCCTGAAGTCTCTGATTAAGTACAAAGGCTATCAGGTGGCTCCCGCTGAATTGGAATCCATCTTGCTCCAACACCCCAACATCTTCGA',
'CGCAGGTGTCGCAGGTCTTCCCGACGATGACGCCGGTGAACTTCCCGCCGCCGTTGTTGTTTTGGAGCACGGAAAGACGATGACGGAAAAAGAG',
'ATCGTGGATTACGTCGCCAGTCAAGTAACAACCGCGAAAAAGTTGCGCGGAGGAGTTGTGTTTGTGGACGAAGTACCGAAAGGTCTTACCGGAA',
'AACTCGACGCAAGAAAAATCAGAGAGATCCTCATAAAGGCCAAGAAGGGCGGAAAGATCGCCGTG')

## Execute ModCon
finalSequence <- ModCon(cds, 1001)

## Print final cds sequence with the alternative SD nucleotide surrounding
print(finalSequence)




## More parameters can be set for use of ModCon when not optimizing to 100\% (e.g. 50\%)

## Execute ModCon
finalSequence <- ModCon(cds, 1001, upChangeCodonsIn=16, downChangeCodonsIn=16,
modconMode='optimalContext', sdMaximalHBS=10, acMaximalMaxent=4,
optiRate=50, nGenerations=5, parentSize=200, startParentSize=800,
bestRate=50, semiLuckyRate=10, luckyRate=5, mutationRate=1e-03, nCores=1)

## Print final cds sequence with the alternative SD nucleotide surrounding
print(finalSequence)

}
