# Function code

increaseGTsiteStrength <- function(cds, sdSeqStartPosition) {

  ## Find position of last nucleotid in first codon covering SD
  startPos <- sdSeqStartPosition
  while (startPos%%3 != 0) startPos <- startPos - 1
  startPos <- startPos + 1

  ## Find position of last nucleotid in last codon covering SD
  endPos <- sdSeqStartPosition + 10
  while (endPos%%3 != 0) endPos <- endPos + 1

  startCodonNr <- (startPos + 2)/3
  endCodonNr <- (endPos)/3

  functionSeq <- substr(cds, startPos, endPos)

  ## Create Codon Matrix by splitting up the sequence by 3nt
  fanFunc <- createCodonMatrix(functionSeq)

  ## Define wt and HEXCO sequence
  after <- paste(fanFunc[1, ], collapse="")

  ## Get position of gt or gc in either of the sequences
  gtPos <- gregexpr(pattern="GT", after)[[1]]
  gtPos <- gtPos[(gtPos <= nchar(after) - 7) & (gtPos >= 4)]
  gtPosSeq <- substr(rep(after, length(gtPos)), gtPos - 3, gtPos + 7)
  gtPosSeq <- hbg$hbs[match(gtPosSeq, hbg$seq)]

  sdExchange <- data.frame(gtPosi=gtPos, hbs=gtPosSeq)

  if (nrow(sdExchange) > 0) {

    ## Get the number of the affected surrounding codons
    sdExchange$upstreamCodonToChange <- ceiling(sdExchange$gtPosi/3) - 1
    sdExchange$downstreamCodonToChange <- sdExchange$upstreamCodonToChange +
      3

    ## Get the respective Codon Code

    sdExchange$codon_code <- paste(Codons$AA[match(fanFunc[2,
                                  (sdExchange$upstreamCodonToChange[1]):
                                    (sdExchange$downstreamCodonToChange[1])],
                                                   Codons$seq)], collapse="-")


    sdExchange$old_seq <- after

    sdExchange <- na.omit(sdExchange)


    ## Create list with position and related dataframe for exchange called 'suche'
    sdList <- split(sdExchange, seq(nrow(sdExchange)))
    sdList <- lapply(sdList, function(x) {
      AA <- strsplit(as.character(x$codon_code), "-")[[1]]
      insertOpt <- expand.grid(c(Codons$seq[Codons$AA == AA[1]]),
                               c(Codons$seq[Codons$AA == AA[2]]),
                               c(Codons$seq[Codons$AA == AA[3]]),
                               c(Codons$seq[Codons$AA == AA[4]]))

      insertOpt$code_seq <- paste0(insertOpt[, 1], insertOpt[, 2],
                                   insertOpt[, 3], insertOpt[, 4])
      outlist <- data.frame(code_seq=insertOpt$code_seq)
      outlist$SD_ID <- x$SD_ID
      outlist$SD_hbs <- x$hbs
      outlist$sum_hex <- as.numeric(lapply(as.character(outlist$code_seq),
                                           calculateHZEIint))
      outlist$downstreamCodonToChange <- x$downstreamCodonToChange
      outlist$upstreamCodonToChange <- x$upstreamCodonToChange
      outlist$pos <- x$gtPosi
      outlist$old_seq <- x$old_seq

      outlist
    })

    ## Get new HBS seq
    sdList <- data.frame(rbindlist(sdList))
    sdList$hbs_seq_search <- as.character(sdList$code_seq)

    ## Get maximal HBS
    sdList$hbs_a <- lapply(sdList$hbs_seq_search, function(x) {

      positions <- gregexpr(pattern="GT", x)[[1]]
      positions <- positions[(positions > 3) & (positions <= nchar(x) - 7)]
      positionsSeq <- substr(rep(x, length(positions)), positions - 3, positions +
                               7)
      positionsSeq <- hbg$hbs[match(positionsSeq, hbg$seq)]

      if (length(positionsSeq) != 0)
        maxHBS <- max(positionsSeq) else maxHBS <- 0
      if (length(positionsSeq) != 0)
        meanHBS <- mean(positionsSeq) else meanHBS <- 0
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
    sdList$sum_hex_abs <- abs(sdList$sum_hex)

    sdList <- sdList[order(sdList$n_hbs, sdList$max_hbs, -sdList$sum_hex_abs,
                           decreasing=TRUE), ]
    sdList <- sdList[1, ]
    finalInserts <- data.frame(seq=sdList$code_seq, up_cod=sdList$upstreamCodonToChange,
                               down_cod=sdList$downstreamCodonToChange)

    finalInserts$seq <- as.character(finalInserts$seq)

    i <- 1

    fanFunc[2, (finalInserts[i, 2]:
                  finalInserts[i, 3])] <- c(substr(finalInserts[i, 1], 1, 3),
                                            substr(finalInserts[i, 1], 4, 6),
                                            substr(finalInserts[i, 1], 7, 9),
                                            substr(finalInserts[i, 1], 10, 12))


  }


  sdIncreaseInputSeq <- paste(fanFunc[2, ], collapse="")
  cdsSplit <- strsplit(as.character(cds), "")[[1]]
  cdsSplit[startPos:endPos] <- strsplit(as.character(sdIncreaseInputSeq), "")[[1]]
  cdsReturn <- paste(cdsSplit, collapse="")

  return(cdsReturn)



}

