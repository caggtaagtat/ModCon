# Function to turn down 5' splice sites
degradeSDs <- function(fanFunc, maxhbs=10, increaseHZEI=TRUE) {

  ## Check if fanFunc matrix has the right format
  if (!is.logical(increaseHZEI)) {
    stop("ERROR: Setting of variable 'increaseHZEI' not as expected.",
                " Must be one logical value.")
  }

  ## Check if fanFunc matrix has the right format
  if ((dim(fanFunc)[[1]] != 2) | (dim(fanFunc)[[2]] < 4)) {
    stop("ERROR: Dimensions of codon matrix 'fanFunc' not as expected.",
                " Must show two rows and at least 4 columns")
  }

  ## Check if fanFunc matrix has the right format
  if (!is.numeric(maxhbs)) {
    stop("ERROR: Setting of variable 'maxhbs' not as expected. Must",
                " hold one numeric value.")
  }

  ## Check if codonsInsert table is present
  if (!any(ls() == "codonsInsert"))
    codonsInsert <- Codons

  ## Define wt and HEXCO sequence
  after <- paste(fanFunc[1, ], collapse="")

  ## Get position of gt or gc in either of the sequences
  gtPos <- gregexpr(pattern="GT", after)[[1]]
  gtPos <- gtPos[(gtPos <= nchar(after) - 7) & (gtPos >= 4)]
  gtPosSeq <- substr(rep(after, length(gtPos)), gtPos - 3, gtPos + 7)
  gtPosSeq <- hbg$hbs[match(gtPosSeq, hbg$seq)]

    sdExchange <- data.frame(gtPosi=gtPos, hbs=gtPosSeq)
  sdExchange <- sdExchange[sdExchange$hbs > maxhbs, ]

  if (nrow(sdExchange) > 0) {

    ## Get the number of the affected surrounding codons
    sdExchange$upstreamCodonToChange <- ceiling(sdExchange$gtPosi/3) - 1
    sdExchange$downstreamCodonToChange <- sdExchange$upstreamCodonToChange +
      3

    ## Get the respective Codon Code
    sdExchange$codon_code <- lapply(seq_len(nrow(sdExchange)), function(x) {
      paste(Codons$AA[match(fanFunc[2, (sdExchange$upstreamCodonToChange[x]):
                                      (sdExchange$downstreamCodonToChange[x])],
                            Codons$seq)], collapse="-")
    })

    sdExchange$SD_ID <- seq_len(nrow(sdExchange))

    ## Get surrounding sequences
    sdExchange$up_cod <- ""
    sdExchange$up_cod[sdExchange$upstreamCodonToChange > 1] <- fanFunc[1,
                                                                       sdExchange$upstreamCodonToChange[sdExchange$upstreamCodonToChange >
                                                                                                             1] - 1]
    sdExchange$down_cod <- as.character(lapply(sdExchange$downstreamCodonToChange,
                                               function(x) {

                                                 if (x < ((nchar(after)/3) - 3)) {

                                                   paste(fanFunc[1, (x + 1):(x + 3)], collapse="")

                                                 } else {

                                                   if (x != (nchar(after)/3)) paste(fanFunc[1, (x + 1):
                                                                                    (nchar(after)/3)],
                                                                                    collapse="") else ""

                                                 }


                                               }))

    sdExchange$old_seq <- as.character(lapply(seq_len(nrow(sdExchange)), function(x) {

      paste(fanFunc[1, (sdExchange$upstreamCodonToChange[x]:
                          sdExchange$downstreamCodonToChange[x])],  collapse="")

    }))

    sdExchange <- na.omit(sdExchange)


    ## Create list with position and related dataframe for exchange called 'suche'
    sdList <- split(sdExchange, seq(nrow(sdExchange)))
    sdList <- lapply(sdList, function(x) {
      AA <- strsplit(as.character(x$codon_code), "-")[[1]]
      insertOpt <- expand.grid(c(codonsInsert$seq[codonsInsert$AA == AA[1]]),
                               c(codonsInsert$seq[codonsInsert$AA == AA[2]]),
                               c(codonsInsert$seq[codonsInsert$AA == AA[3]]),
                               c(codonsInsert$seq[codonsInsert$AA == AA[4]]))
      insertOpt$code_seq <- paste0(insertOpt[, 1],
                                   insertOpt[, 2],
                                   insertOpt[, 3],
                                   insertOpt[, 4])
      outlist <- data.frame(code_seq=insertOpt$code_seq)
      outlist$sum_hex <- as.numeric(lapply(as.character(outlist$code_seq), calculateHZEIint))

      ## If overall HZEI should be decreased, multiply calculated HZEIint by -1
      if (!increaseHZEI)
        outlist$sum_hex <- outlist$sum_hex * -1
      outlist$SD_ID <- x$SD_ID
      outlist$SD_hbs <- x$hbs
      outlist$downstreamCodonToChange <- x$downstreamCodonToChange
      outlist$upstreamCodonToChange <- x$upstreamCodonToChange
      outlist$pos <- x$gtPosi
      outlist$up_cod <- x$up_cod
      outlist$down_cod <- x$down_cod
      outlist$old_seq <- x$old_seq

      outlist
    })

    ## Get new HBS seq
    sdList <- data.frame(rbindlist(sdList))
    sdList$hbs_seq_search <- paste0(sdList$up_cod, sdList$code_seq, sdList$down_cod)

    ## Get maximal HBS
    sdList$hbs_a <- lapply(sdList$hbs_seq_search, function(x) {

      positions <- gregexpr(pattern="GT", x)[[1]]
      positions <- positions[(positions > 3) & (positions <= nchar(x) - 7)]
      positionSeq <- substr(rep(x, length(positions)), positions - 3, positions +
                              7)
      positionSeq <- hbg$hbs[match(positionSeq, hbg$seq)]

      if (length(positionSeq) != 0)
        maxHBS <- max(positionSeq) else maxHBS <- 0
      if (length(positionSeq) != 0)
        meanHBS <- mean(positionSeq) else meanHBS <- 0
      nHBS <- length(positions)

      paste(maxHBS, nHBS, meanHBS)
    })

    sdList$hbs_a <- as.character(sdList$hbs_a)

    ## Get max HBS
    sdList$max_hbs <- as.numeric(lapply(seq_len(nrow(sdList)), function(x) {
      strsplit(as.character(sdList$hbs_a[x]), " ")[[1]][[1]]
    }))

    ## Get number of SDs
    sdList$n_hbs <- as.numeric(lapply(seq_len(nrow(sdList)), function(x) {
      strsplit(as.character(sdList$hbs_a[x]), " ")[[1]][[2]]
    }))

    ## Get mean HBS
    sdList$mean_hbs <- as.numeric(lapply(seq_len(nrow(sdList)), function(x) {
      strsplit(as.character(sdList$hbs_a[x]), " ")[[1]][[3]]
    }))
    sdList$hbs_a <- NULL

    ## Get backup surrounding
    sdList$old <- 0
    sdList$old[sdList$code_seq == sdList$old_seq] <- 1

    sdList$old_search_seq <- paste0(sdList$up_cod, sdList$old_seq, sdList$down_cod)

    sdList$old_seq_HZEI <- as.numeric(lapply(sdList$old_search_seq, calculateHZEIint))
    sdList$new_seq_HZEI <- as.numeric(lapply(sdList$hbs_seq_search, calculateHZEIint))
    sdList$HZEI_diff <- sdList$new_seq_HZEI - sdList$old_seq_HZEI

    ## If overall HZEI should be decreased, multiply calculated HZEI_diff by -1
    if (!increaseHZEI)
      sdList$HZEI_diff <- sdList$HZEI_diff * -1

    sdList$max_hbs_group <- cut(sdList$max_hbs, breaks=c(-1, 3, 6, 9, 12, 15,
                                                           18, 21, 24))
    sdList$mean_hbs_group <- cut(sdList$mean_hbs, breaks=c(-1, 3, 6, 9, 12,
                                                             15, 18, 21, 24))

    ## Exclude some replacement candidates based on HBS and HZEI
    finalInserts <- lapply(unique(sdList$SD_ID), function(x) {

      subList <- sdList[sdList$SD_ID == x, ]
      subList <- subList[order(subList$max_hbs_group, subList$n_hbs, subList$mean_hbs_group,
                               -subList$HZEI_diff, decreasing=FALSE), ]
      subList <- subList[1, ]

      data.frame(seq=subList$code_seq, up_cod=subList$upstreamCodonToChange,
                 down_cod=subList$downstreamCodonToChange)

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
  }

  return(paste(fanFunc[2, ], collapse=""))

}


