# ModCon
ModCon is a tool to manipulate the nucleotide sequence surrounding of splice donors, in order to either enhance or silence their usage. It is based on predictions of the HEXplorer tool, which estimates the property of a genomic sequences, to recruit splice regulatory proteins.
<br/><br/>
## Installation and Usage
#### 1. Install dependencies
The ModCon tool was developed as an R package (version ≥ R-4.0.3) using integrated Perl scripts (version ≥ Perl 5.18.1) for MaxEntScan score calculations. Therefore, if necessary please install Perl and R first. If wanted, the script performs parallel HZEI calculation parallel with multiple cores. In total three R packages are required to execute the application, namely “parallel”, “utils” and "data.table".
<br/><br/>
In order to execute the ModCon R-script, please download every file of the ModCon repository for installation or use the devtools option to install R packages from github.<br/>
`library(devtools)`<br/>
`install_github("caggtaagtat/ModCon")`
<br/><br/>
#### 2. Define coding sequence
First, please define the coding sequence of interest which holds the splice donor in question. The first nucleotide of the entered sequence must be the first nucleotide of a codon, in order to correctly read the encoded amino acid sequence.<br/>
`cds <- "AGTAACAGTTAGACCAAAGGATAGA..." # coding sequence holding the SD of interest`
<br/><br/>
#### 3. Enter position of first SD sequence nucleotide
Now, state the position of the first nucleotide of the splice donor sequence within the entered coding sequence, again be redefining the exemplary variable sdSeqStart.<br/>
`sdSeqStart <- 402`
<br/><br/>
#### 4. Execute the main ModCon function
Select whether to optimize SD context for SD usage enhancement or silencing. Per default, the script will generate the optimal (or suboptimal) sequence surrounding for SD usage enhancement (or silencing).<br/>
`finalSequence <- ModCon(cds, sdSeqStart)` or <br/>
`finalSequence <- ModCon(cds, sdSeqStart, optimizeContext=FALSE)`
<br/><br/>
#### 5. OPTIONAL: Adjust additional parameters
Additionally, if the circumstances may require adjustment of further parameters, please do so.<br/>
`finalSequence <- ModCon(cds, sdSeqStartPosition, upChangeCodonsIn=16, downChangeCodonsIn=16, optimizeContext="optimalContext", sdMaximalHBS=10, acMaximalMaxent=4, optiRate=100, nGenerations=30, parentSize=300, startParentSize=1000, bestRate=40, semiLuckyRate=20, luckyRate=5, mutationRate=1e-04, nCores=-1)`
<br/><br/>

## Parameter settings
| Parameters             | Default           | Definition  |
|:----------------------:|:-----------------:| :-----------------------------------------------------------------------------|
|cds                     | CDS of interest   | Coding sequence which holds the splice donor in question.             |
|sdSeqStartPosition      | 1001              | Position of the first nucleotide of the splice donor sequence within the entered codong sequence in column CDS                       |
|upChangeCodonsIn        | 16                | Number of codons to be substituted upstream of the splice donor sequence.            |
|downChangeCodonsIn      | 16                | Number of codons to be substituted downstream of the splice donor sequence.              |
|optimizeContext         | TRUE              | Define, whether the script should either adjust the SD surrounding by codon selection, in order to enhance (TRUE, default) or silence (FALSE) splice donor usage.  |
|sdMaximalHBS            | 10                | Splice donors within the output sequence, which show a HBS of X or higher, will be attempted to be degraded by codon selection.                              |
|acMaximalMaxent         | 4                 | Splice acceptors within the output sequence, which show a MaxEntScan score of X or higher, will be attempted to be degraded by codon selection.             |
|optiRate                | 100               | Change the extent the splice site HEXplorer weight of the donor site should be adjusted. |
|nGenerations            | 50                | Maximal numbers of generations to generate, before the loop of sequence recombination is terminated, in case the condition to end the loop will never be reached. |
|parentSize              | 300               | Number of sequences in the parental generation, which comes from sequences of either the initial or the filial sequence population.    |
|startParentSize         | 1000              | Number of sequences which are randomly generated during the generation of the initial parent population.                                  |
|bestRate                | 40                | X % of parental sequences with the highest fitness are selected to make up parts of the reproducing population.                      |
|semiLuckyRate           | 20                | X % of parental sequences are selected to make up parts of the reproducing population, based on a probability which is derived from the sequence fitness.  |
|luckyRate               | 5                 | X % of sequences are randomly selected to make up parts of the reproducing population |
|mutationRate            | 1e-04             | With a chance of X % each codon within a filial sequence can mutate to a different codon, encoding the same amino acid.                |
