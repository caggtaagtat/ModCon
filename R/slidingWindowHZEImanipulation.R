## Creating sequence population
slidingWindowHZEImanipulation <- function(inSeq,
                                          increaseHZEI=TRUE) {

  ## Check if codonsInsert table is present
  if (!any(ls() == "codonsInsert"))
    codonsInsert <- Codons

  ## increaseHZEI
  if (!is.logical(increaseHZEI)) {
    stop("ERROR: Setting of variable 'increaseHZEI' not as",
                " expected. Must be one logical value.")
  }

  if (!all(strsplit(inSeq, "")[[1]] %in% c("a", "c", "g", "t", "G", "C", "T", "A")))
    stop("ERROR during setting of variable 'CDS'. The entered sequence must be a",
         " character string of A C G and T.")

  ## Transformation of the sequence to upper case
  wtSequence <- toupper(inSeq)

  ## Create Codon Matrix by splitting up the sequence by 3nt
  fan <- createCodonMatrix(wtSequence)
  fan[2, ] <- Codons$AA[match(fan[1, ], Codons$seq)]

  ## The first step is special and will be taken outside the following loop
  start <- 1
  end <- 4

  ## Calculate every possibility to encode the first 4 amino acids
  AA <- fan[2, c(start:end)]
  newSequences <- expand.grid(c(codonsInsert$seq[codonsInsert$AA == AA[1]]),
                              c(codonsInsert$seq[codonsInsert$AA == AA[2]]),
                              c(codonsInsert$seq[codonsInsert$AA == AA[3]]),
                              c(codonsInsert$seq[codonsInsert$AA == AA[4]]))
  newSequences$seq <- apply(newSequences, 1, paste, collapse="")

  ## Calculate HZEI integral for every sequnence variant
  newSequences$HZEI_int <- as.numeric(lapply(newSequences$seq,
                                             calculateHZEIint))

  ## Access last hexamer of each sequence, since hexamers consititute a stable
  ## boarder in of the overlapping effect of HZEI integral calculation
  newSequences$last_hexamer <- substr(newSequences$seq,
                                      (nchar(newSequences$seq) - 5),
                                      nchar(newSequences$seq))

  ## Get maximal or minimal HZEI integral per last hexamer
  if (increaseHZEI)
    bestHZEIperHexamer <- tapply(newSequences$HZEI_int,
                                 newSequences$last_hexamer, max)
  if (!increaseHZEI)
    bestHZEIperHexamer <- tapply(newSequences$HZEI_int,
                                 newSequences$last_hexamer, min)

  ## For every sequence pool with the same last hexamer, select the one with the
  ## highest HZEI integral
  newSequences$max_gr <- bestHZEIperHexamer[match(newSequences$last_hexamer,
                                                  names(bestHZEIperHexamer))]
  newSequences$max <- 0
  newSequences$max[newSequences$HZEI_int == newSequences$max_gr] <- 1
  newSequences <- newSequences[newSequences$max == 1, ]

  newSequences$Var2 <- NULL
  newSequences$Var3 <- NULL
  newSequences$Var4 <- NULL

  start <- end + 1
  finalSequence <- ""

  ## Repeat the same steps as above for the next 4 amino acids
  for (i in 2:floor(length(fan[2, ])/4)) {

    ## Go to the next 4 amino acids and save every potential codon combination
    end <- i * 4
    AA <- fan[2, c(start:end)]
    newSequences <- expand.grid(newSequences$seq,
                                c(codonsInsert$seq[codonsInsert$AA == AA[1]]),
                                c(codonsInsert$seq[codonsInsert$AA == AA[2]]),
                                c(codonsInsert$seq[codonsInsert$AA == AA[3]]),
                                c(codonsInsert$seq[codonsInsert$AA == AA[4]]))

    newSequences$seq <- apply(newSequences, 1, paste, collapse="")
    newSequences$Var1 <- NULL
    newSequences$Var2 <- NULL
    newSequences$Var3 <- NULL
    newSequences$Var4 <- NULL
    newSequences$Var5 <- NULL
    newSequences$HZEI_int <- as.numeric(lapply(newSequences$seq, calculateHZEIint))
    newSequences$last_hexamer <- substr(newSequences$seq, (nchar(newSequences$seq) -
                                                             5), nchar(newSequences$seq))

    if (increaseHZEI)
      bestHZEIperHexamer <- tapply(newSequences$HZEI_int, newSequences$last_hexamer,
                                   max)
    if (!increaseHZEI)
      bestHZEIperHexamer <- tapply(newSequences$HZEI_int, newSequences$last_hexamer,
                                   min)

    newSequences$max_gr <- bestHZEIperHexamer[match(newSequences$last_hexamer,
                                                    names(bestHZEIperHexamer))]
    newSequences$max <- 0
    newSequences$max[newSequences$HZEI_int == newSequences$max_gr] <- 1
    newSequences <- newSequences[newSequences$max == 1, ]

    if (end == length(fan[2, ])) {

      if (increaseHZEI)
        newSequences <- newSequences[order(newSequences$max_gr, decreasing=TRUE),
                                     ]
      if (!increaseHZEI)
        newSequences <- newSequences[order(newSequences$max_gr, decreasing=FALSE),
                                     ]
      finalSequence <- paste0(finalSequence, newSequences$seq[1])

    } else {

      if (length(unique(substr(newSequences$seq, 1, (nchar(newSequences$seq) -
                                                     12)))) == 1) {
        finalSequence <- paste0(finalSequence,
                                unique(substr(newSequences$seq, 1,
                                              (nchar(newSequences$seq) - 12))))
        newSequences$seq <- substr(newSequences$seq,
                                   (nchar(newSequences$seq) - 11),
                                   (nchar(newSequences$seq)))
      }
    }

    start <- end + 1
  }

  ## In case the length of the amino acids sequence is not divisible by 4 do the
  ## same for the remaining amino acids
  if (end < length(fan[2, ])) {

    end <- length(fan[2, ])

    AA <- fan[2, c(start:end)]
    newSequences <- expand.grid(newSequences$seq,
                                c(codonsInsert$seq[codonsInsert$AA == AA[1]]),
                                c(codonsInsert$seq[codonsInsert$AA == AA[2]]),
                                c(codonsInsert$seq[codonsInsert$AA == AA[3]]))
    newSequences$seq <- apply(newSequences, 1, paste, collapse="")
    newSequences$Var1 <- NULL
    newSequences$Var2 <- NULL
    newSequences$Var3 <- NULL
    newSequences$Var4 <- NULL
    newSequences$seq <- gsub("NA", "", newSequences$seq)
    newSequences <- unique(newSequences)
    newSequences$HZEI_int <- as.numeric(lapply(newSequences$seq, calculateHZEIint))
    newSequences$last_hexamer <- substr(newSequences$seq,
                                        (nchar(newSequences$seq) - 5),
                                        nchar(newSequences$seq))

    if (increaseHZEI)
      bestHZEIperHexamer <- tapply(newSequences$HZEI_int,
                                   newSequences$last_hexamer, max)
    if (!increaseHZEI)
      bestHZEIperHexamer <- tapply(newSequences$HZEI_int,
                                   newSequences$last_hexamer, min)

    newSequences$max_gr <- bestHZEIperHexamer[match(newSequences$last_hexamer,
                                                    names(bestHZEIperHexamer))]
    newSequences$max <- 0
    newSequences$max[newSequences$HZEI_int == newSequences$max_gr] <- 1
    newSequences <- newSequences[newSequences$max == 1, ]

    if (increaseHZEI)
      newSequences <- newSequences[order(newSequences$HZEI_int, decreasing=TRUE),
                                   ]
    if (!increaseHZEI)
      newSequences <- newSequences[order(newSequences$HZEI_int, decreasing=FALSE),
                                   ]



    finalSequence <- paste0(finalSequence, newSequences$seq[1])
  }

  ## Save the final sequence
  return(finalSequence)

}
