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
#### 2. Enter coding sequence
To replace the example sequence with the coding sequence of interest which holds the splice donor in question, redefine the internal variable CDS. The first nucleotide of the entered sequence must be the first nucleotide of a codon, in order to correctly read the encoded amino acid sequence.<br/>
`CDS <- "AGTAACAGTTAGACCAAAGGATAGA..." # coding sequence holding the SD of interest`
<br/><br/>
#### 3. Enter position of first SD sequence nucleotide
After replacing the exemplary coding sequence within the script, please state the position of the first nucleotide of the splice donor sequence within the entered coding sequence, again be redefining the exemplary variable SD_seq_start_position.<br/>
`SD_seq_start_position <- 402`
<br/><br/>
#### 4. Select program type
Select whether to optimize SD context for SD usage enhancement or silencing. Per default, the script will generate the optimal sequence surrounding for SD usage enhancement.<br/>
`Programm_version <- "Optimal_context"` or <br/>
`Programm_version <- "Sub-optimal_context"`
<br/><br/>
#### 5. Adjust parameters
Additionally, if the circumstances may require adjustment of further parameters, please do so carefully. 
<br/><br/>
#### 6. Load dependencies
Than execute the ModCon function. The variables up_change_codons and down_change_codons, specifying the length of sequence adjustment are set to a default of 16.<br/>
`final_sequence <- ModCon(CDS, SD_seq_start_position, up_change_codons, down_change_codons, Programm_version )`
<br/><br/>

## Parameter settings
                                               
| Parameters             | Default           | Definition  |
|:----------------------:|:-----------------:| :-----------------------------------------------------------------------------|
|n_generations           | 50                | Maximal numbers of generations to generate, before the loop of sequence recombination is terminated, in case the condition to end the loop will never be reached. |
| parent_size            | 300               |   Number of sequences in the parental generation, which comes from sequences of either the initial or the filial sequence population.    |
| START_parent_size      | 1000              |  Number of sequences which are randomly generated during the generation of the initial parent population.                                  |
|Best_rate               | 40                | X % of parental sequences with the highest fitness are selected to make up parts of the reproducing population.                      |
| Semi_Lucky_rate        | 20                | X % of parental sequences are selected to make up parts of the reproducing population, based on a probability which is derived from the sequence fitness.  |
| Lucky_rate             | 5                 | X % of sequences are randomly selected to make up parts of the reproducing population |
|Mutation_rate           | 1e-05             | With a chance of X % each codon within a filial sequence can mutate to a different codon, encoding the same amino acid.                |
| opti_rate              | 100               | Since hundreds of sequences are generated with increasingly optimized properties, it is possible to set the desired degree of maximization or minimization in percent, with 100% representing the highest or lowest possible. |
| SDMaximalHBS           | 10                | Splice donors within the output sequence, which show a HBS of X or higher, will be attempted to be degraded.                              |
|ACMaximalMaxent         | 4                 | Splice acceptors within the output sequence, which show a MaxEntScan score of X or higher, will be attempted to be degraded.            |
|up_change_codons         | 16                 | Number of codons to be substituted upstream of the splice donor sequence.            |
|down_change_codons         | 16                 | Number of codons to be substituted downstream of the splice donor sequence.              |
| CDS                    | CDS of Luciferase | Coding sequence which holds the splice donor in question.             |
|SD_seq_start_position   | 1001              | Position of the first nucleotide of the splice donor sequence within the entered codong sequence in column CDS                       |
|Programm_version        | Optimal_context   | Define, whether the script should either adjust the SD surrounding in order to enhance or silence splice donor usage. Alternative input would be "Sub-optimal_context".               |
