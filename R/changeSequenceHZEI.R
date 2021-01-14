## Creating sequence population
changeSequenceHZEI <- function(inSeq, increaseHZEI=TRUE, nGenerations=30,
                               parentSize=300, startParentSize=1000,
                               bestRate=50, semiLuckyRate=20, luckyRate=5,
                               mutationRate=1e-04, optiRate=100, sdMaximalHBS=10,
                               acMaximalMaxent=4, nCores=-1) {

  ## Error messages for wrong entries sdMaximalHBS
  if (sdMaximalHBS < 1.8 | sdMaximalHBS > 23.8) {
    stop(paste0("ERROR during setting of variable 'sdMaximalHBS'.",
                " HBS limit must be within 1.8 to 23.8."))
  }

  ## sdMaximalHBS
  if (!is.numeric(sdMaximalHBS)) {
    stop(paste0("ERROR during setting of variable 'acMaximalMaxent'.",
                " MaxEntScan score limit must be numeric."))
  }

  ## optiRate
  if (!is.numeric(optiRate) | optiRate < 0) {
    stop(paste0("ERROR during setting of variable 'optiRate'.",
                " optiRate must be numeric and not lower than 0 ."))
  }

  ## nGenerations
  if (!is.numeric(nGenerations) | nGenerations < 2) {
    stop(paste0("ERROR during setting of variable 'nGenerations'.",
                " nGenerations must be numeric and greater 1."))
  }

  ## parentSize
  if (!is.numeric(parentSize) | parentSize < 10) {
    stop(paste0("ERROR during setting of variable 'parentSize'.",
                " parentSize must be numeric and greater 9."))
  }

  ## startParentSize
  if (!is.numeric(startParentSize) | startParentSize < 10) {
    stop(paste0("ERROR during setting of variable 'startParentSize'.",
                " startParentSize must be numeric and greater 9."))
  }

  ## bestRate
  if (!is.numeric(bestRate) | bestRate <= 0) {
    stop(paste0("ERROR during setting of variable 'bestRate'.",
                " bestRate must be numeric and greater 0."))
  }

  ## semiLuckyRate
  if (!is.numeric(semiLuckyRate) | semiLuckyRate <= 0) {
    stop(paste0("ERROR during setting of variable 'semiLuckyRate'.",
                " semiLuckyRate must be numeric and greater 0."))
  }

  ## luckyRate
  if (!is.numeric(luckyRate) | luckyRate <= 0) {
    stop(paste0("ERROR during setting of variable 'luckyRate'.",
                " luckyRate must be numeric and greater 0."))
  }

  ## mutationRate
  if (!is.numeric(mutationRate) | mutationRate < 0 | mutationRate > 1) {
    stop(paste0("ERROR during setting of variable 'mutationRate'.",
                " mutationRate must be numeric and range from 0 to 1."))
  }

  ## increaseHZEI
  if (!is.logical(increaseHZEI)) {
    stop(paste0("ERROR: Setting of variable 'increaseHZEI' not as",
                " expected. Must be one logical value."))
  }

  ## nCores
  if (!is.numeric(nCores)) {
    stop(paste0("ERROR: Setting of variable 'nCores' not as expected.",
                " Must be one numeric value e.g. 1 or 4."))
  }

  wtSequence <- inSeq
  ## Check for correct syntax
  if (!all(strsplit(wtSequence, "")[[1]] %in% c("a", "c", "g", "t", "G", "C",
                                                "T", "A")))
    stop("The entered sequence must be a character string of A C G and T.")


  ## Transformation to upper case
  wildtypeSequenz <- toupper(wtSequence)


  ## If optimization rate is not 100, apply the genetic algorithm for a heuristic
  ## calculation of the sequence with the entered optimization rate
  if (optiRate != 100) {

    ## Create Codon Matrix by splitting up the sequence by 3nt
    fan <- createCodonMatrix(wildtypeSequenz)
    fan[2, ] <- Codons$AA[match(fan[1, ], Codons$seq)]

    ## Setup cluster for parralell computation with as many cores as possible
    if(nCores == -1){
      nCoresIntern <- detectCores()
    }else{
      nCoresIntern <- nCores
    }
    clust <- makeCluster(nCoresIntern)

    ## Setting itiration number
    numberGenerations <- nGenerations

    ## Define 3 parent populations
    parentsize <- parentSize
    startingPopulation <- startParentSize

    ## Create codonsInsert object, to allow differences between lookup table Codons
    ## and insertion table codonsInsert
    if (!any(ls() == "codonsInsert"))
      codonsInsert <- Codons

    ## Export data and functions to clusters
    parent0 <- character()
    clusterExport(clust, list("generateRandomCodonsPerAA", "fan"), envir=environment())
    clusterExport(clust, list("getOverlappingVectorsFromVector", "calculateHZEIint",
                              "hex", "Codons", "startingPopulation", "parent0",
                              "increaseHZEI", "hex2", "wtSequence", "codonsInsert"),
                  envir=environment())

    ## Generate first parent population
    parent0 <- as.character(parLapply(clust, seq_len(startingPopulation), function(x) {
      parent0[[x]] <- generateRandomCodonsPerAA(fan[2, ])
    }))

    ## Do pre-optimization with sliding window
    seqSpikes <- slidingWindowHZEImanipulation(wildtypeSequenz, increaseHZEI=increaseHZEI)
    parent0[(startingPopulation + 1):(startingPopulation + (startingPopulation *
                                                              0.1))] <- seqSpikes

    ## Genetic Algorithm Convert 3 different values for population generation to
    ## percentages
    whoMatesBestPercent <- bestRate
    whoMatesSemiRandom <- ceiling((semiLuckyRate/100) * parentsize)
    whoMatesLuckily <- ceiling((luckyRate/100) * parentsize)
    mutationChance <- mutationRate

    ## Generate 3 output lists
    best <- numeric()
    means <- numeric()
    saveAllSequences <- character()

    ## Define generation iteration variable
    generation <- 1
    STATUS <- "RUNNING"

    ## As long as STATUS implies, continue to generate new generations
    while (STATUS == "RUNNING") {
      
      ## If all parent sequences are the same stop the Genetic algorithm
      if(length(unique(parent0)) == 1) STATUS <- "STOP"

      ## Select fittest sequence plus some random individuals for recombination
      matingIndividuals <- selectMatingIndividuals(parent0,
                                                   whoMatesBestPercent=whoMatesBestPercent,
                                                   whoMatesSemiRandom=whoMatesSemiRandom,
                                                   whoMatesLuckily=whoMatesLuckily,
                                                   clust, increaseHZEI=increaseHZEI)

      ## Create next generation by recombination
      childrens <- createFilialSequencePopulation(matingIndividuals, parentsize)

      ## Mutate the filial sequence
      childrens <- mutatePopulation(childrens, mutationChance)

      ## Save all sequences
      saveAllSequences <- c(saveAllSequences, parent0)

      ## Use the new generation as the new parent generation for the next iteration
      parent0 <- childrens

      ## Save best and mean
      status <- selectBestAndMean(parent0, clust, increaseHZEI=increaseHZEI)
      best <- c(best, status[1])
      means <- c(means, status[2])

      ## Proceed to next generation by adding one to the iterating variable generation
      generation <- generation + 1

      ## If generation greater 30 or maximal number, check if last 20 generation showed
      ## no significant improvement of best value and if so end the loop
      if (generation > 30) {
        bestVar <- best[(length(best) - 20):length(best)]
        bestVar <- bestVar - best[length(best)]

        if ((generation > 30) & (all(bestVar < 2 & bestVar > -2)))
          STATUS <- "STOP"

      }

      if (generation >= numberGenerations)
        STATUS <- "STOP"

      ## Report the current generation number and best HZEI integral
      message(paste("Best HZEI-integral in generation ", generation, ": ",
                    status[1]))

    }

    ## Create result list
    results <- list()

    ## Save highes achieved HZEI-score integral
    results[[1]] <- max(best)
    if (!increaseHZEI)
      results[[1]] <- min(best)

    ## Adding the wildtype sequence and calculating the HZEI-score integrals
    saveAllSequences <- c(saveAllSequences, wildtypeSequenz)

    ## Generate Dataframe
    candidates <- data.frame(seq=saveAllSequences, HZEIint=0)

    ## make the table unique
    candidates <- unique(candidates)

    ## Calculate integral
    candidates$HZEIint <- as.numeric(parLapply(clust, as.character(candidates$seq),
                                               calculateHZEIint))

    ## Stop the Cluster
    stopCluster(clust)

    ## Calculate optimization rank and HZEI-integral change in percentage
    candidates$type <- "new"
    candidates$type[nrow(candidates)] <- "wt"
    candidates$HZEIint <- as.numeric(candidates$HZEIint)
    candidates <- candidates[order(candidates$HZEIint, decreasing=TRUE), ]
    if (!increaseHZEI)
      candidates <- candidates[order(candidates$HZEIint, decreasing=FALSE), ]
    candidates$type[1] <- "max"
    if (!increaseHZEI)
      candidates$type[1] <- "min"
    candidates$HZEI_change <- candidates$HZEIint - candidates$HZEIint[candidates$type ==
                                                                        "wt"]
    candidates$HZEI_change <- round((candidates$HZEI_change/candidates$HZEI_change[1]) *
                                      100, 4)

    results[[2]] <- candidates

    ## Select customly optimized sequence
    selectedSeq <- candidates
    selectedSeq <- selectedSeq[selectedSeq$HZEI_change <= optiRate, ]
    selectedSeq <- selectedSeq[order(selectedSeq$HZEI_change, decreasing=TRUE),
                               ]
    selectedSeq <- as.character(selectedSeq$seq[[1]])

    ## If optimization rate is 100, apply the sliding window algorithm
  } else {

    ## Apply the fast and exact window optimization method
    selectedSeq <- slidingWindowHZEImanipulation(wildtypeSequenz,
                                                 increaseHZEI=increaseHZEI)

    ## Create result list
    results <- list()

    ## Save highes achieved HZEI-score integral
    results[[1]] <- selectedSeq

    ## Adding the wildtype sequence and calculating the HZEI-score integrals
    saveAllSequences <- c(selectedSeq, wildtypeSequenz)

    ## Generate Dataframe of unique seqs
    candidates <- data.frame(seq=saveAllSequences, HZEIint=0)
    candidates <- unique(candidates)

    ## Calculate integral per seq
    candidates$HZEIint <- lapply(as.character(candidates$seq),
                                 calculateHZEIint)

    ## Calculate optimization rank and change in percent
    candidates$type <- "new"
    candidates$type[candidates$seq == wildtypeSequenz] <- "wt"
    candidates$HZEIint <- as.numeric(candidates$HZEIint)
    candidates <- candidates[order(candidates$HZEIint, decreasing=TRUE), ]
    if (!increaseHZEI)
      candidates <- candidates[order(candidates$HZEIint, decreasing=FALSE), ]
    candidates$type[1] <- "max"
    if (!increaseHZEI)
      candidates$type[1] <- "min"
    candidates$HZEI_change <- candidates$HZEIint - candidates$HZEIint[candidates$type ==
                                                                        "wt"]
    candidates$HZEI_change <- round((candidates$HZEI_change/candidates$HZEI_change[1]) *
                                      100, 4)
    results[[2]] <- candidates

  }

  ## Degrade SDs with HBS greater than the allowed maximal SD strength Select
  ## sequence
  optSeq <- selectedSeq

  ## Create Codon Matrix by splitting up the sequence by 3nt
  fan <- createCodonMatrix(optSeq)

  ## Degrade SDs
  optSeq <- degradeSDs(fan, maxhbs=sdMaximalHBS, increaseHZEI=increaseHZEI)
  optSeq <- as.character(optSeq)

  ## Create Codon Matrix by splitting up the sequence by 3nt
  fan <- createCodonMatrix(optSeq)

  ## Degrade SAs keeping SDs in mind
  optSeqFinal <- degradeSAs(fan, maxhbs=sdMaximalHBS, maxME=acMaximalMaxent,
                            increaseHZEI=increaseHZEI)
  results[[3]] <- optSeqFinal

  ## Convert seqs to character values
  candidates$seq <- as.character(candidates$seq)

  ## Calculate HZEI intergral change in percent
  candidates$HZEI_change <- round(candidates$HZEI_change/10)
  ct <- data.frame(tapply(candidates$HZEIint, candidates$HZEI_change, max))
  if (!increaseHZEI)
    ct <- data.frame(tapply(candidates$HZEIint, candidates$HZEI_change, min))
  ct$opti_level <- as.numeric(row.names(ct)) * 10
  ct <- ct[order(ct$opti_level, decreasing=TRUE), ]
  names(ct) <- c("HZEI_int", "opti_level")
  row.names(ct) <- NULL

  ## Save the sequences
  ct$sequence <- ""
  ct$sequence <- candidates$seq[match(as.numeric(ct$HZEI_int),
                                      as.numeric(candidates$HZEIint))]
  ct$sequence[ct$opti_level == 0] <- wildtypeSequenz
  ct$HZEI_int[ct$opti_level == 0] <- calculateHZEIint(as.character(wildtypeSequenz))
  results[[4]] <- ct

  ## Return the results list
  return(results)

}


