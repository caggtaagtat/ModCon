## Recombination of mating sequences
createFilialSequencePopulation <- function(sequenceVector,
                                           generateNrecombinedSequences){

  filialPopulation <- character()

  nPerPair <- ceiling(generateNrecombinedSequences/(length(sequenceVector)/2))

  for (i in seq_len(length(sequenceVector)/2)) {
    for (j in seq_len(nPerPair)) {
      filialPopulation <- c(filialPopulation,
                            recombineTwoSequences(sequenceVector[i],
                            sequenceVector[(length(sequenceVector) - i)],
                            c(0.1, 0.6, 0.3)))

    }

  }

  filialPopulation <- filialPopulation[seq_len(generateNrecombinedSequences)]

  return(filialPopulation)

}
