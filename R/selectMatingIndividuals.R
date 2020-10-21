# Select fittest plus some random individuals for recombination
selectMatingIndividuals <- function(inputGeneration,
                                    whoMatesBestPercent=40,
                                    whoMatesSemiRandom=20,
                                    whoMatesLuckily=5,
                                    clust,
                                    increaseHZEI=TRUE) {

  inputGeneration <- data.frame(seq=as.character(inputGeneration), fitt=0)
  inputGeneration$fitt <- as.numeric(parLapply(clust,
                                               as.character(inputGeneration$seq),
                                               calculateHZEIint))
  if (!increaseHZEI)
    inputGeneration$fitt <- inputGeneration$fitt * -1
  inputGeneration <- inputGeneration[order(inputGeneration$fitt,
                                           decreasing=TRUE),
                                     ]
  inputGenerationMates <- head(inputGeneration,
                               (whoMatesBestPercent/100) * nrow(inputGeneration))
  inputGenerationMates2 <- inputGeneration[sample(seq_len(length(inputGeneration)),
                                                  whoMatesLuckily, replace=TRUE), ]
  inputGenerationMates <- rbind(inputGenerationMates, inputGenerationMates2)
  inputGenerationMates <- inputGenerationMates$seq
  inputGeneration$fitt <- inputGeneration$fitt + abs(min(inputGeneration$fitt))
  inputGenerationMates3 <- sample(inputGeneration$seq, whoMatesSemiRandom,
                                  prob=inputGeneration$fitt)
  inputGenerationMates <- c(as.character(inputGenerationMates),
                            as.character(inputGenerationMates3))
  return(inputGenerationMates)
}
