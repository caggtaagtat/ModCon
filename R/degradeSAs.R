## Function to degrade strength of 3' splice sites
degradeSAs <- function(fanFunc, maxhbs = 10, maxME = 4, increaseHZEI = TRUE) {

  ## Check if fanFunc matrix has the right format
  if ((dim(fanFunc)[[1]] != 2) | (dim(fanFunc)[[2]] < 4)) {
    stop(paste0("ERROR: Dimensions of codon matrix 'fanFunc' not as expected.",
                " Must show two rows and at least 4 columns"))
  }

  ## maxhbs
  if (!is.numeric(maxhbs)) {
    stop(paste0("ERROR: Setting of variable 'maxhbs' not as expected.",
                " Must hold one numeric value."))
  }

  ## increaseHZEI
  if (!is.logical(increaseHZEI)) {
    stop(paste0("ERROR: Setting of variable 'increaseHZEI' not as expected.",
                " Must be one logical value."))
  }

  ## maxME
  if (!is.numeric(maxME)) {
    stop(paste0("ERROR: Setting of variable 'maxME' not as expected. Must ",
                "hold one numeric value."))
  }

  ## Check if codonsInsert table is present
  if (!any(ls() == "codonsInsert"))
    codonsInsert <- Codons

  ## Define wt and HEXCO sequence
  after <- paste(fanFunc[1, ], collapse = "")

  ## Finding acceptor possitions with a maxentscore greater zero
  agPos <- gregexpr(pattern = "AG", after)[[1]]
  agPos <- agPos[(agPos <= nchar(after) - 5) & (agPos > 18)]


  ## In case there is at least one SA to be degraded
  if (length(agPos) > 0) {
    agSeq <- as.character(lapply(agPos, function(x) {
      substr(after, x - 18, x + 4)
    }))

    agSeqMax <- calculateMaxEntScanScore(agSeq, 3)

    ## Create dataframe with SA positions which have to be degraded
    saExchange <- data.frame(agPosi = agPos, ME = as.numeric(agSeqMax))
    saExchange <- saExchange[saExchange$ME > maxME, ]

    ## In case there is at least one SA to be degraded
    if (nrow(saExchange) > 0) {

      ## Get the number of the affected surrounding codons
      saExchange$upstreamCodonToChange <- ceiling(saExchange$agPosi/3) - 2
      saExchange$downstreamCodonToChange <- saExchange$upstreamCodonToChange +
        3

      ## Get the respective Codon Code
      saExchange$codon_code <- lapply(seq_len(nrow(saExchange)), function(x) {
        paste(Codons$AA[match(fanFunc[2, (saExchange$upstreamCodonToChange[x]):
                                        (saExchange$downstreamCodonToChange[x])],
                              Codons$seq)], collapse = "-")
      })

      saExchange$SD_ID <- seq_len(nrow(saExchange))

      ## Get surrounding sequences
      saExchange$up_cod <- ""
      saExchange$up_cod[saExchange$upstreamCodonToChange > 1] <- fanFunc[1,
                                                                  saExchange$upstreamCodonToChange[saExchange$upstreamCodonToChange > 1] - 1]

      saExchange$down_cod <- as.character(lapply(saExchange$downstreamCodonToChange,
                                                 function(x) {

                                                   if (x < ((nchar(after)/3) - 3)) {

                                                     paste(fanFunc[1, (x + 1):(x + 3)], collapse = "")

                                                   } else {

                                                     if (x != (nchar(after)/3)) paste(fanFunc[1, (x + 1):
                                                                                          (nchar(after)/3)],
                                                                                      collapse = "") else ""

                                                   }

                                                 }))


      ## Save old sequence, in case degradation is not properly possible
      saExchange$old_seq <- as.character(lapply(seq_len(nrow(saExchange)), function(x) {

        paste(fanFunc[1, (saExchange$upstreamCodonToChange[x]:
                            saExchange$downstreamCodonToChange[x])], collapse = "")

      }))


      ## Create list with position and related dataframe for exchange called 'suche'
      saList <- split(saExchange, seq(nrow(saExchange)))
      saList <- lapply(saList, function(x) {
        AA <- strsplit(as.character(x$codon_code), "-")[[1]]
        insertOpt <- expand.grid(c(codonsInsert$seq[codonsInsert$AA == AA[1]]),
                                 c(codonsInsert$seq[codonsInsert$AA == AA[2]]),
                                 c(codonsInsert$seq[codonsInsert$AA == AA[3]]),
                                 c(codonsInsert$seq[codonsInsert$AA == AA[4]]))

        insertOpt$code_seq <- paste0(insertOpt[, 1],
                                     insertOpt[, 2],
                                     insertOpt[, 3],
                                     insertOpt[, 4])
        outlist <- data.frame(code_seq = insertOpt$code_seq)
        outlist$code_seq <- as.character(outlist$code_seq)
        outlist$sum_hex <- as.numeric(lapply(outlist$code_seq, calculateHZEIint))
        ## If overall HZEI should be decreased, multiply calculated HZEIint by -1
        if (!increaseHZEI)
          outlist$sum_hex <- outlist$sum_hex * -1
        outlist$SD_ID <- x$SD_ID
        outlist$SD_hbs <- x$hbs
        outlist$downstreamCodonToChange <- x$downstreamCodonToChange
        outlist$upstreamCodonToChange <- x$upstreamCodonToChange
        outlist$pos <- x$agPosi
        outlist$up_cod <- x$up_cod
        outlist$down_cod <- x$down_cod
        outlist$old_seq <- x$old_seq

        outlist
      })

      ## Get new HBS seq
      saList <- data.frame(rbindlist(saList))
      saList$hbs_seq_search <- paste0(saList$up_cod, saList$code_seq, saList$down_cod)

      ## Get maximal HBS
      saList$hbs_a <- lapply(saList$hbs_seq_search, function(x) {

        positions <- gregexpr(pattern = "GT", x)[[1]]
        positions <- positions[(positions > 3) & (positions <= nchar(x) - 7)]
        positions_seq <- substr(rep(x, length(positions)), positions - 3, positions +
                                  7)
        positions_seq <- hbg$hbs[match(positions_seq, hbg$seq)]

        if (length(positions_seq) != 0)
          maxHBS <- max(positions_seq) else maxHBS <- 0
        if (length(positions_seq) != 0)
          meanHBS <- mean(positions_seq) else meanHBS <- 0
        nHBS <- length(positions)

        paste(maxHBS, nHBS, meanHBS)
      })

      saList$hbs_a <- as.character(saList$hbs_a)

      ## Get max HBS
      saList$max_hbs <- as.numeric(lapply(seq_len(nrow(saList)), function(x) {
        strsplit(saList$hbs_a[x], " ")[[1]][[1]]
      }))

      ## Get number of SDs
      saList$n_hbs <- as.numeric(lapply(seq_len(nrow(saList)), function(x) {
        strsplit(saList$hbs_a[x], " ")[[1]][[2]]
      }))

      ## Get mean HBS
      saList$mean_hbs <- as.numeric(lapply(seq_len(nrow(saList)), function(x) {
        strsplit(saList$hbs_a[x], " ")[[1]][[3]]
      }))
      saList$hbs_a <- NULL

      ## Get backup surrounding
      saList$old <- 0
      saList$old[saList$code_seq == saList$old_seq] <- 1

      saList$old_search_seq <- paste0(saList$up_cod, saList$old_seq, saList$down_cod)

      saList$old_seq_HZEI <- as.numeric(lapply(saList$old_search_seq, calculateHZEIint))
      saList$new_seq_HZEI <- as.numeric(lapply(saList$hbs_seq_search, calculateHZEIint))

      saList$HZEI_diff <- saList$new_seq_HZEI - saList$old_seq_HZEI

      ## If overall HZEI should be decreased, multiply calculated HZEI_diff by -1
      if (!increaseHZEI)
        saList$HZEI_diff <- saList$HZEI_diff * -1

      ## Group HBS values per exchange
      saList$max_hbs_group <- cut(saList$max_hbs, breaks = c(-1, 3, 6, 9, 12, 15,
                                                             18, 21, 24))
      saList$mean_hbs_group <- cut(saList$mean_hbs, breaks = c(-1, 3, 6, 9, 12,
                                                               15, 18, 21, 24))

      ## Find AG dinucleotides
      saList$old_AGs <- as.numeric(lapply(saList$old_search_seq, function(x) {
        length(as.numeric(gregexpr("AG", x)[[1]]))
      }))
      saList$new_AGs <- as.numeric(lapply(saList$hbs_seq_search, function(x) {
        length(as.numeric(gregexpr("AG", x)[[1]]))
      }))
      saList$diff_AGs <- saList$new_AGs - saList$old_AGs

      ## Exclude some replacement candidates based on HBS and HZEI
      finalInserts <- lapply(unique(saList$SD_ID), function(x) {

        subList <- saList[saList$SD_ID == x, ]
        subList <- subList[order(subList$diff_AGs, subList$max_hbs_group, subList$n_hbs,
                                 subList$mean_hbs_group, -subList$HZEI_diff, decreasing = FALSE), ]
        subList <- subList[1, ]

        data.frame(seq = subList$code_seq, up_cod = subList$upstreamCodonToChange,
                   down_cod = subList$downstreamCodonToChange)

      })

      finalInserts <- data.frame(rbindlist(finalInserts))
      finalInserts$seq <- as.character(finalInserts$seq)

      ## Replace the codons to replace
      for (i in seq_len(nrow(finalInserts))) {

        fanFunc[2, (finalInserts[i, 2]:
                      finalInserts[i, 3])] <- c(substr(finalInserts[i, 1], 1, 3),
                                                substr(finalInserts[i, 1], 4, 6),
                                                substr(finalInserts[i, 1], 7, 9),
                                                substr(finalInserts[i, 1], 10, 12))

      }
    }}

  ## Return the final sequence as a character string
  return(paste(fanFunc[2, ], collapse = ""))

}

